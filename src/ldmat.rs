//! LD matrix storage. Full row-major f64 (reference), and packed upper-triangle (symmetric) storage in f64 or
//! f32: each stored element is read once per product and used twice (y_i += a_ij x_j, y_j += a_ij x_i), so
//! memory traffic is halved; f32 halves it again. Accumulation is always f64.
use crate::simd;

pub trait LdMat: Sync {
    fn p(&self) -> usize;
    fn get(&self, i: usize, j: usize) -> f64;
    /// out = R v
    fn matvec(&self, v: &[f64], out: &mut [f64]);
}

pub struct LdF64 { pub p: usize, pub r: Vec<f64> }
impl LdMat for LdF64 {
    fn p(&self) -> usize { self.p }
    #[inline] fn get(&self, i: usize, j: usize) -> f64 { self.r[i * self.p + j] }
    fn matvec(&self, v: &[f64], out: &mut [f64]) {
        let p = self.p;
        for (i, o) in out.iter_mut().enumerate() { *o = simd::dot(&self.r[i * p..(i + 1) * p], v) }
    }
}

/// Packed upper triangle, row i holds R[i, i..p].
pub struct LdSym<T> { pub p: usize, pub r: Vec<T>, off: Vec<usize> }
fn offsets(p: usize) -> Vec<usize> { let mut o = Vec::with_capacity(p + 1); let mut s = 0; for i in 0..p { o.push(s); s += p - i } o.push(s); o }

impl LdMat for LdSym<f64> {
    fn p(&self) -> usize { self.p }
    #[inline] fn get(&self, i: usize, j: usize) -> f64 { let (i, j) = if i <= j { (i, j) } else { (j, i) }; self.r[self.off[i] + j - i] }
    fn matvec(&self, v: &[f64], out: &mut [f64]) {
        out.iter_mut().for_each(|o| *o = 0.0);
        for i in 0..self.p {
            let row = &self.r[self.off[i]..self.off[i + 1]];
            // diagonal + upper part: out[i] += row . v[i..];  out[i+1..] += row[1..] * v[i]
            out[i] += simd::dot(row, &v[i..]);
            simd::axpy(v[i], &row[1..], &mut out[i + 1..]);
        }
    }
}
impl LdMat for LdSym<f32> {
    fn p(&self) -> usize { self.p }
    #[inline] fn get(&self, i: usize, j: usize) -> f64 { let (i, j) = if i <= j { (i, j) } else { (j, i) }; self.r[self.off[i] + j - i] as f64 }
    fn matvec(&self, v: &[f64], out: &mut [f64]) {
        out.iter_mut().for_each(|o| *o = 0.0);
        for i in 0..self.p {
            let row = &self.r[self.off[i]..self.off[i + 1]];
            out[i] += simd::dot_f32(row, &v[i..]);
            simd::axpy_f32(v[i], &row[1..], &mut out[i + 1..]);
        }
    }
}

fn read_raw(path: &str, p: usize) -> Result<Vec<u8>, String> {
    let b = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    if b.len() != p * p * 8 { return Err(format!("{path}: {} bytes, expected {} for p = {p}", b.len(), p * p * 8)) }
    Ok(b)
}
#[inline] fn at(b: &[u8], k: usize) -> f64 { let x = f64::from_le_bytes(b[k * 8..k * 8 + 8].try_into().unwrap()); if x.is_finite() { x } else { 0.0 } }

pub fn read_full(path: &str, p: usize) -> Result<LdF64, String> {
    let b = read_raw(path, p)?;
    Ok(LdF64 { p, r: (0..p * p).map(|k| at(&b, k)).collect() })
}
/// Symmetrised packed copy: (R + R')/2 on the upper triangle, as susieR forces symmetry.
pub fn read_sym_f64(path: &str, p: usize) -> Result<LdSym<f64>, String> {
    let b = read_raw(path, p)?; let off = offsets(p);
    let mut r = Vec::with_capacity(off[p]);
    for i in 0..p { for j in i..p { r.push(if i == j { at(&b, i * p + j) } else { (at(&b, i * p + j) + at(&b, j * p + i)) / 2.0 }) } }
    Ok(LdSym { p, r, off })
}
pub fn read_sym_f32(path: &str, p: usize) -> Result<LdSym<f32>, String> {
    let s = read_sym_f64(path, p)?;
    Ok(LdSym { p, r: s.r.iter().map(|x| *x as f32).collect(), off: s.off })
}

/// Low-rank factor from reference genotypes: R = W W' with W = X_std / sqrt(N) (p x N, row-major).
/// R v = W (W' v) costs 2pN instead of p^2, and W stays in cache for typical reference panels (N ~ 500).
pub struct LdFactor<T> { pub p: usize, pub n: usize, pub w: Vec<T> }
impl LdMat for LdFactor<f64> {
    fn p(&self) -> usize { self.p }
    fn get(&self, i: usize, j: usize) -> f64 { let n = self.n; simd::dot(&self.w[i * n..(i + 1) * n], &self.w[j * n..(j + 1) * n]) }
    fn matvec(&self, v: &[f64], out: &mut [f64]) {
        let n = self.n; let mut t = vec![0.0f64; n];
        for i in 0..self.p { if v[i] != 0.0 { simd::axpy(v[i], &self.w[i * n..(i + 1) * n], &mut t) } }
        for (i, o) in out.iter_mut().enumerate() { *o = simd::dot(&self.w[i * n..(i + 1) * n], &t) }
    }
}
impl LdMat for LdFactor<f32> {
    fn p(&self) -> usize { self.p }
    fn get(&self, i: usize, j: usize) -> f64 { let n = self.n; self.w[i * n..(i + 1) * n].iter().zip(&self.w[j * n..(j + 1) * n]).map(|(a, b)| *a as f64 * *b as f64).sum() }
    fn matvec(&self, v: &[f64], out: &mut [f64]) {
        let n = self.n; let mut t = vec![0.0f64; n];
        for i in 0..self.p { if v[i] != 0.0 { simd::axpy_f32(v[i], &self.w[i * n..(i + 1) * n], &mut t) } }
        for (i, o) in out.iter_mut().enumerate() { *o = simd::dot_f32(&self.w[i * n..(i + 1) * n], &t) }
    }
}
/// geno.bin: p*N little-endian f64 of standardised genotypes (rows = variants), N given.
pub fn read_factor(path: &str, p: usize, n: usize) -> Result<LdFactor<f64>, String> {
    let b = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    if b.len() != p * n * 8 { return Err(format!("{path}: {} bytes, expected {}", b.len(), p * n * 8)) }
    let s = 1.0 / (n as f64).sqrt();
    Ok(LdFactor { p, n, w: b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap()) * s).collect() })
}
pub fn read_factor_f32(path: &str, p: usize, n: usize) -> Result<LdFactor<f32>, String> {
    let f = read_factor(path, p, n)?; Ok(LdFactor { p, n, w: f.w.iter().map(|x| *x as f32).collect() })
}

/// Packed symmetric R built from a genotype factor (for small p, where p^2/2 < 2pN).
pub fn sym_from_factor(f: &LdFactor<f64>) -> LdSym<f64> {
    use rayon::prelude::*;
    let (p, n) = (f.p, f.n); let off = offsets(p);
    let rows: Vec<Vec<f64>> = (0..p).into_par_iter().map(|i| (i..p).map(|j| simd::dot(&f.w[i * n..(i + 1) * n], &f.w[j * n..(j + 1) * n])).collect()).collect();
    LdSym { p, r: rows.concat(), off }
}
