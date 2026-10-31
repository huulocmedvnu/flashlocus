//! Dot-product kernels with runtime CPU dispatch (AVX-512F, AVX2+FMA, scalar). One binary runs everywhere;
//! on banhbo (Xeon 8168, Skylake-SP) the AVX-512 path is used.
use std::sync::OnceLock;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

static LEVEL: OnceLock<u8> = OnceLock::new();
pub fn level() -> u8 {
    *LEVEL.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        {
            if std::env::var("FLASHLOCUS_NO_AVX512").is_err() && is_x86_feature_detected!("avx512f") { return 2 }
            if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") { return 1 }
        }
        0
    })
}
pub fn level_name() -> &'static str { match level() { 2 => "avx512f", 1 => "avx2+fma", _ => "scalar" } }

#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    #[cfg(target_arch = "x86_64")]
    unsafe {
        match level() { 2 => return dot_avx512(a, b), 1 => return dot_avx2(a, b), _ => {} }
    }
    dot_scalar(a, b)
}

/// f32 row (LD stored in single precision) times f64 vector, accumulated in f64.
#[inline]
pub fn dot_f32(a: &[f32], b: &[f64]) -> f64 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        match level() { 2 => return dot_f32_avx512(a, b), 1 => return dot_f32_avx2(a, b), _ => {} }
    }
    a.iter().zip(b).map(|(x, y)| *x as f64 * y).sum()
}

/// y += a * x
#[inline]
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    unsafe { match level() { 2 => return axpy_avx512(a, x, y), 1 => return axpy_avx2(a, x, y), _ => {} } }
    for (yi, xi) in y.iter_mut().zip(x) { *yi += a * xi }
}
#[inline]
pub fn axpy_f32(a: f64, x: &[f32], y: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    unsafe { if level() == 2 { return axpy_f32_avx512(a, x, y) } }
    for (yi, xi) in y.iter_mut().zip(x) { *yi += a * *xi as f64 }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn axpy_avx512(a: f64, x: &[f64], y: &mut [f64]) {
    let n = x.len(); let (px, py) = (x.as_ptr(), y.as_mut_ptr()); let va = _mm512_set1_pd(a);
    let mut i = 0;
    while i + 16 <= n {
        _mm512_storeu_pd(py.add(i), _mm512_fmadd_pd(va, _mm512_loadu_pd(px.add(i)), _mm512_loadu_pd(py.add(i))));
        _mm512_storeu_pd(py.add(i + 8), _mm512_fmadd_pd(va, _mm512_loadu_pd(px.add(i + 8)), _mm512_loadu_pd(py.add(i + 8))));
        i += 16;
    }
    while i + 8 <= n { _mm512_storeu_pd(py.add(i), _mm512_fmadd_pd(va, _mm512_loadu_pd(px.add(i)), _mm512_loadu_pd(py.add(i)))); i += 8 }
    if i < n {
        let k: __mmask8 = ((1u16 << (n - i)) - 1) as __mmask8;
        _mm512_mask_storeu_pd(py.add(i), k, _mm512_fmadd_pd(va, _mm512_maskz_loadu_pd(k, px.add(i)), _mm512_maskz_loadu_pd(k, py.add(i))));
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn axpy_avx2(a: f64, x: &[f64], y: &mut [f64]) {
    let n = x.len(); let (px, py) = (x.as_ptr(), y.as_mut_ptr()); let va = _mm256_set1_pd(a);
    let mut i = 0;
    while i + 4 <= n { _mm256_storeu_pd(py.add(i), _mm256_fmadd_pd(va, _mm256_loadu_pd(px.add(i)), _mm256_loadu_pd(py.add(i)))); i += 4 }
    while i < n { *py.add(i) += a * *px.add(i); i += 1 }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn axpy_f32_avx512(a: f64, x: &[f32], y: &mut [f64]) {
    let n = x.len(); let (px, py) = (x.as_ptr(), y.as_mut_ptr()); let va = _mm512_set1_pd(a);
    let mut i = 0;
    while i + 8 <= n { _mm512_storeu_pd(py.add(i), _mm512_fmadd_pd(va, _mm512_cvtps_pd(_mm256_loadu_ps(px.add(i))), _mm512_loadu_pd(py.add(i)))); i += 8 }
    while i < n { *py.add(i) += a * *px.add(i) as f64; i += 1 }
}

/// (max, sum_j exp(x_j - max)): log-sum-exp pieces with a vectorised exp (AVX-512: range reduction to
/// 2^k * exp(r), |r| <= ln2/2, degree-12 Taylor polynomial, relative error ~1e-16).
pub fn lse_parts(x: &[f64]) -> (f64, f64) {
    let m = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() { return (m, 0.0) }
    #[cfg(target_arch = "x86_64")]
    unsafe { if level() == 2 { return (m, sum_exp_shift_avx512(x, m)) } }
    (m, x.iter().map(|v| (v - m).exp()).sum())
}
/// out_j = exp(x_j - m)
pub fn exp_shift(x: &[f64], m: f64, out: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    unsafe { if level() == 2 { return exp_shift_avx512(x, m, out) } }
    for (o, v) in out.iter_mut().zip(x) { *o = (v - m).exp() }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn exp_pd(x: __m512d) -> __m512d {
    let x = _mm512_max_pd(x, _mm512_set1_pd(-745.0));
    let k = _mm512_roundscale_pd::<0x08>(_mm512_mul_pd(x, _mm512_set1_pd(std::f64::consts::LOG2_E)));
    let r = _mm512_fnmadd_pd(k, _mm512_set1_pd(6.93147180369123816490e-01), x);
    let r = _mm512_fnmadd_pd(k, _mm512_set1_pd(1.90821492927058770002e-10), r);
    const C: [f64; 13] = [1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0, 1.0 / 5040.0, 1.0 / 40320.0,
                          1.0 / 362880.0, 1.0 / 3628800.0, 1.0 / 39916800.0, 1.0 / 479001600.0];
    let mut p = _mm512_set1_pd(C[12]);
    for i in (0..12).rev() { p = _mm512_fmadd_pd(p, r, _mm512_set1_pd(C[i])) }
    _mm512_scalef_pd(p, k)
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn sum_exp_shift_avx512(x: &[f64], m: f64) -> f64 {
    let n = x.len(); let px = x.as_ptr(); let vm = _mm512_set1_pd(m);
    let (mut a0, mut a1) = (_mm512_setzero_pd(), _mm512_setzero_pd());
    let mut i = 0;
    while i + 16 <= n {
        a0 = _mm512_add_pd(a0, exp_pd(_mm512_sub_pd(_mm512_loadu_pd(px.add(i)), vm)));
        a1 = _mm512_add_pd(a1, exp_pd(_mm512_sub_pd(_mm512_loadu_pd(px.add(i + 8)), vm)));
        i += 16;
    }
    let mut s = _mm512_reduce_add_pd(_mm512_add_pd(a0, a1));
    while i < n { s += (*px.add(i) - m).exp(); i += 1 }
    s
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn exp_shift_avx512(x: &[f64], m: f64, out: &mut [f64]) {
    let n = x.len(); let (px, po) = (x.as_ptr(), out.as_mut_ptr()); let vm = _mm512_set1_pd(m);
    let mut i = 0;
    while i + 8 <= n { _mm512_storeu_pd(po.add(i), exp_pd(_mm512_sub_pd(_mm512_loadu_pd(px.add(i)), vm))); i += 8 }
    while i < n { *po.add(i) = (*px.add(i) - m).exp(); i += 1 }
}

fn dot_scalar(a: &[f64], b: &[f64]) -> f64 {
    let (mut s0, mut s1, mut s2, mut s3) = (0.0, 0.0, 0.0, 0.0);
    let m = a.len() / 4 * 4;
    let mut i = 0;
    while i < m { s0 += a[i] * b[i]; s1 += a[i + 1] * b[i + 1]; s2 += a[i + 2] * b[i + 2]; s3 += a[i + 3] * b[i + 3]; i += 4 }
    let mut s = (s0 + s1) + (s2 + s3);
    while i < a.len() { s += a[i] * b[i]; i += 1 }
    s
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn dot_avx512(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len(); let (pa, pb) = (a.as_ptr(), b.as_ptr());
    let (mut c0, mut c1, mut c2, mut c3) = (_mm512_setzero_pd(), _mm512_setzero_pd(), _mm512_setzero_pd(), _mm512_setzero_pd());
    let mut i = 0;
    while i + 32 <= n {
        c0 = _mm512_fmadd_pd(_mm512_loadu_pd(pa.add(i)), _mm512_loadu_pd(pb.add(i)), c0);
        c1 = _mm512_fmadd_pd(_mm512_loadu_pd(pa.add(i + 8)), _mm512_loadu_pd(pb.add(i + 8)), c1);
        c2 = _mm512_fmadd_pd(_mm512_loadu_pd(pa.add(i + 16)), _mm512_loadu_pd(pb.add(i + 16)), c2);
        c3 = _mm512_fmadd_pd(_mm512_loadu_pd(pa.add(i + 24)), _mm512_loadu_pd(pb.add(i + 24)), c3);
        i += 32;
    }
    while i + 8 <= n { c0 = _mm512_fmadd_pd(_mm512_loadu_pd(pa.add(i)), _mm512_loadu_pd(pb.add(i)), c0); i += 8 }
    if i < n {
        let k: __mmask8 = ((1u16 << (n - i)) - 1) as __mmask8;
        c1 = _mm512_fmadd_pd(_mm512_maskz_loadu_pd(k, pa.add(i)), _mm512_maskz_loadu_pd(k, pb.add(i)), c1);
    }
    _mm512_reduce_add_pd(_mm512_add_pd(_mm512_add_pd(c0, c1), _mm512_add_pd(c2, c3)))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn dot_avx2(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len(); let (pa, pb) = (a.as_ptr(), b.as_ptr());
    let (mut c0, mut c1, mut c2, mut c3) = (_mm256_setzero_pd(), _mm256_setzero_pd(), _mm256_setzero_pd(), _mm256_setzero_pd());
    let mut i = 0;
    while i + 16 <= n {
        c0 = _mm256_fmadd_pd(_mm256_loadu_pd(pa.add(i)), _mm256_loadu_pd(pb.add(i)), c0);
        c1 = _mm256_fmadd_pd(_mm256_loadu_pd(pa.add(i + 4)), _mm256_loadu_pd(pb.add(i + 4)), c1);
        c2 = _mm256_fmadd_pd(_mm256_loadu_pd(pa.add(i + 8)), _mm256_loadu_pd(pb.add(i + 8)), c2);
        c3 = _mm256_fmadd_pd(_mm256_loadu_pd(pa.add(i + 12)), _mm256_loadu_pd(pb.add(i + 12)), c3);
        i += 16;
    }
    let s = _mm256_add_pd(_mm256_add_pd(c0, c1), _mm256_add_pd(c2, c3));
    let mut t = [0.0f64; 4]; _mm256_storeu_pd(t.as_mut_ptr(), s);
    let mut r = (t[0] + t[1]) + (t[2] + t[3]);
    while i < n { r += a[i] * b[i]; i += 1 }
    r
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn dot_f32_avx512(a: &[f32], b: &[f64]) -> f64 {
    let n = a.len(); let (pa, pb) = (a.as_ptr(), b.as_ptr());
    let (mut c0, mut c1) = (_mm512_setzero_pd(), _mm512_setzero_pd());
    let mut i = 0;
    while i + 16 <= n {
        c0 = _mm512_fmadd_pd(_mm512_cvtps_pd(_mm256_loadu_ps(pa.add(i))), _mm512_loadu_pd(pb.add(i)), c0);
        c1 = _mm512_fmadd_pd(_mm512_cvtps_pd(_mm256_loadu_ps(pa.add(i + 8))), _mm512_loadu_pd(pb.add(i + 8)), c1);
        i += 16;
    }
    let mut r = _mm512_reduce_add_pd(_mm512_add_pd(c0, c1));
    while i < n { r += *pa.add(i) as f64 * *pb.add(i); i += 1 }
    r
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn dot_f32_avx2(a: &[f32], b: &[f64]) -> f64 {
    let n = a.len(); let (pa, pb) = (a.as_ptr(), b.as_ptr());
    let (mut c0, mut c1) = (_mm256_setzero_pd(), _mm256_setzero_pd());
    let mut i = 0;
    while i + 8 <= n {
        c0 = _mm256_fmadd_pd(_mm256_cvtps_pd(_mm_loadu_ps(pa.add(i))), _mm256_loadu_pd(pb.add(i)), c0);
        c1 = _mm256_fmadd_pd(_mm256_cvtps_pd(_mm_loadu_ps(pa.add(i + 4))), _mm256_loadu_pd(pb.add(i + 4)), c1);
        i += 8;
    }
    let s = _mm256_add_pd(c0, c1); let mut t = [0.0f64; 4]; _mm256_storeu_pd(t.as_mut_ptr(), s);
    let mut r = (t[0] + t[1]) + (t[2] + t[3]);
    while i < n { r += *pa.add(i) as f64 * *pb.add(i); i += 1 }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn kernels_agree() {
        for n in [1usize, 7, 8, 9, 31, 33, 100, 1001] {
            let a: Vec<f64> = (0..n).map(|i| ((i * 37 % 11) as f64 - 5.0) / 3.0).collect();
            let b: Vec<f64> = (0..n).map(|i| ((i * 13 % 7) as f64 - 3.0) / 7.0).collect();
            let r = dot_scalar(&a, &b);
            assert!((dot(&a, &b) - r).abs() < 1e-10, "n={n}");
            let af: Vec<f32> = a.iter().map(|x| *x as f32).collect();
            assert!((dot_f32(&af, &b) - r).abs() < 1e-4, "f32 n={n}");
            let xs: Vec<f64> = (0..n).map(|i| -700.0 + 1400.0 * i as f64 / n as f64 - (i % 3) as f64 * 0.37).collect();
            let (m, se) = lse_parts(&xs); let se0: f64 = xs.iter().map(|v| (v - m).exp()).sum();
            assert!(((se - se0) / se0).abs() < 1e-14, "lse n={n}");
            let mut o = vec![0.0; n]; exp_shift(&xs, m, &mut o);
            assert!(o.iter().zip(&xs).all(|(u, v)| { let t = (v - m).exp(); t == 0.0 && *u < 1e-300 || ((u - t) / t).abs() < 2e-15 }), "exp n={n}");
            let mut y1 = b.clone(); let mut y2 = b.clone(); axpy(0.7, &a, &mut y1);
            for (u, v) in y2.iter_mut().zip(&a) { *u += 0.7 * v }
            assert!(y1.iter().zip(&y2).all(|(u, v)| (u - v).abs() < 1e-12), "axpy n={n}");
        }
    }
}
