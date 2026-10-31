//! SuSiE-RSS following susieR 0.14.2: susie_rss(z, R, n) -> susie_suff_stat with XtX = (n-1)R,
//! Xty = sqrt(n-1) z_adj, yty = n-1, standardize = TRUE, residual variance fixed, prior variance by
//! optim/Brent on log V in [-30, 15]. The per-effect products XtX b_l are cached, so each IBSS sweep needs
//! L matrix-vector products and the ELBO is O(Lp); the arithmetic is otherwise the same as susieR.
use crate::brent;
use crate::ldmat::LdMat;
use crate::simd::{self, dot};

#[derive(Clone)]
#[allow(dead_code)]
pub struct Params {
    pub l: usize,
    pub scaled_prior_variance: f64,
    pub tol: f64,
    pub max_iter: usize,
    pub coverage: f64,
    pub min_abs_corr: f64,
    pub prior_tol: f64,
    pub check_null_threshold: f64,
    pub estimate_prior_variance: bool,
    pub check_prior: bool,
    pub n_purity: usize,
    pub threads_in_locus: bool,
}
impl Default for Params {
    fn default() -> Self {
        Params { l: 10, scaled_prior_variance: 0.2, tol: 1e-3, max_iter: 100, coverage: 0.95, min_abs_corr: 0.5,
                 prior_tol: 1e-9, check_null_threshold: 0.0, estimate_prior_variance: true, check_prior: true,
                 n_purity: 100, threads_in_locus: false }
    }
}

pub struct CredibleSet { pub effect: usize, pub members: Vec<usize>, pub coverage: f64, pub purity: [f64; 3] }

#[allow(dead_code)] // full susieR fit object kept for library use
pub struct Fit {
    pub p: usize, pub l: usize,
    pub alpha: Vec<f64>, pub mu: Vec<f64>, pub mu2: Vec<f64>, pub lbf_variable: Vec<f64>,
    pub v: Vec<f64>, pub lbf: Vec<f64>, pub kl: Vec<f64>, pub elbo: Vec<f64>,
    pub niter: usize, pub converged: bool, pub pip: Vec<f64>, pub sets: Vec<CredibleSet>,
}

/// Correlation-like matrix in row-major order with a scale: A = s * D^-1 R D^-1.
struct Op<'a, M: LdMat> { r: &'a M, s: f64, inv_csd: Vec<f64> }
impl<'a, M: LdMat> Op<'a, M> {
    fn mul(&self, b: &[f64], out: &mut [f64]) {
        let bb: Vec<f64> = b.iter().zip(&self.inv_csd).map(|(x, c)| x * c).collect();
        self.r.matvec(&bb, out);
        for (o, c) in out.iter_mut().zip(&self.inv_csd) { *o *= self.s * c }
    }
}

fn log_dnorm_diff(bh: f64, v: f64, s2: f64) -> f64 {
    // dnorm(bh, 0, sqrt(V + s2), log) - dnorm(bh, 0, sqrt(s2), log)
    let sd1 = (v + s2).sqrt(); let sd0 = s2.sqrt();
    (-sd1.ln() - 0.5 * (bh / sd1) * (bh / sd1)) - (-sd0.ln() - 0.5 * (bh / sd0) * (bh / sd0))
}

/// Single-effect data. When every shat2 is equal (R with unit diagonal, the usual case) the log BF is
/// a + b * betahat^2 with a, b depending only on V: one fused multiply-add per variant per evaluation.
struct SerData<'a> { bh: Vec<f64>, bh2: Vec<f64>, s2: Vec<f64>, lpw: &'a [f64], s2c: Option<f64>, buf: Vec<f64> }
impl<'a> SerData<'a> {
    fn new(xty: &[f64], d: &[f64], sigma2: f64, lpw: &'a [f64]) -> SerData<'a> {
        let bh: Vec<f64> = xty.iter().zip(d).map(|(x, d)| x / d).collect();
        let s2: Vec<f64> = d.iter().map(|d| sigma2 / d).collect();
        // equal up to rounding (diag of R computed from genotypes is 1 +- 1 ulp): relative tolerance 1e-10
        let s2c = if s2.iter().all(|x| x.is_finite() && ((x - s2[0]) / s2[0]).abs() <= 1e-10) { Some(s2[0]) } else { None };
        let bh2 = bh.iter().map(|x| x * x).collect();
        let p = bh.len();
        SerData { bh, bh2, s2, lpw, s2c, buf: vec![0.0; p] }
    }
    /// fills buf with lpo_j = lbf_j(V) + log prior weight_j
    fn lpo(&mut self, v: f64) {
        if let Some(s2) = self.s2c {
            let sd1 = (v + s2).sqrt(); let sd0 = s2.sqrt();
            let a = -sd1.ln() + sd0.ln();
            let b = 0.5 * (1.0 / (sd0 * sd0) - 1.0 / (sd1 * sd1));
            for j in 0..self.bh.len() { self.buf[j] = a + b * self.bh2[j] + self.lpw[j] }
        } else {
            for j in 0..self.bh.len() {
                self.buf[j] = if self.s2[j].is_infinite() { 0.0 } else { log_dnorm_diff(self.bh[j], v, self.s2[j]) + self.lpw[j] }
            }
        }
    }
    fn loglik(&mut self, v: f64) -> f64 { self.lpo(v); let (m, s) = simd::lse_parts(&self.buf); m + s.ln() }
}

struct Ser { alpha: Vec<f64>, mu: Vec<f64>, mu2: Vec<f64>, lbf: Vec<f64>, v: f64, lbf_model: f64 }

fn ser(xty: &[f64], d: &[f64], v0: f64, sigma2: f64, lpw: &[f64], est: bool, thr: f64) -> Ser {
    let mut sd = SerData::new(xty, d, sigma2, lpw);
    let mut v = v0;
    if est {
        let mut nll = |lv: f64| -sd.loglik(lv.exp());
        let mut lv = brent::fmin(-30.0, 15.0, &mut nll, f64::EPSILON.sqrt());
        if nll(lv) > nll(v.ln()) { lv = v.ln() }
        v = lv.exp();
        if sd.loglik(0.0) + thr >= sd.loglik(v) { v = 0.0 }
    }
    let p = xty.len();
    sd.lpo(v);
    let lpo = std::mem::take(&mut sd.buf);
    let lbf: Vec<f64> = (0..p).map(|j| if sd.s2[j].is_infinite() { 0.0 } else { lpo[j] - lpw[j] }).collect();
    let (mx, sw) = simd::lse_parts(&lpo);
    let mut alpha = vec![0.0; p]; simd::exp_shift(&lpo, mx, &mut alpha);
    for a in alpha.iter_mut() { *a /= sw }
    let mut mu = vec![0.0; p]; let mut mu2 = vec![0.0; p];
    for j in 0..p {
        let pv = 1.0 / (1.0 / v + d[j] / sigma2);
        mu[j] = (1.0 / sigma2) * pv * xty[j];
        mu2[j] = pv + mu[j] * mu[j];
    }
    Ser { alpha, mu, mu2, lbf, v, lbf_model: mx + sw.ln() }
}

/// z: z-scores; r: p x p LD (correlation) row-major; n: GWAS sample size.
pub fn susie_rss<M: LdMat>(z_in: &[f64], r: &M, n: f64, prm: &Params) -> Result<Fit, String> {
    let p = z_in.len();
    if r.p() != p { return Err(format!("LD has p = {}, z has {}", r.p(), p)) }
    if n <= 1.0 { return Err("n must be greater than 1".into()) }
    // susie_rss: PVE adjustment, XtX = (n-1)R, Xty = sqrt(n-1) z, yty = n-1
    let z: Vec<f64> = z_in.iter().map(|x| { let x = if x.is_nan() { 0.0 } else { *x }; ((n - 1.0) / (x * x + n - 2.0)).sqrt() * x }).collect();
    // standardize = TRUE in susie_suff_stat: csd = sqrt(diag(XtX)/(n-1)) = sqrt(R_ii)
    let csd: Vec<f64> = (0..p).map(|i| { let c = r.get(i, i).sqrt(); if c == 0.0 { 1.0 } else { c } }).collect();
    let op = Op { r, s: n - 1.0, inv_csd: csd.iter().map(|c| 1.0 / c).collect() };
    let d: Vec<f64> = (0..p).map(|i| (n - 1.0) * r.get(i, i) / (csd[i] * csd[i])).collect();
    let xty: Vec<f64> = (0..p).map(|i| (n - 1.0).sqrt() * z[i] / csd[i]).collect();
    let yty = n - 1.0;
    let l = prm.l.min(p);
    let sigma2 = yty / (n - 1.0);
    let lpw: Vec<f64> = vec![(1.0 / p as f64 + f64::EPSILON.sqrt()).ln(); p];
    let mut alpha = vec![1.0 / p as f64; l * p];
    let mut mu = vec![0.0; l * p]; let mut mu2 = vec![0.0; l * p];
    let mut lbfv = vec![f64::NAN; l * p];
    let mut v = vec![prm.scaled_prior_variance * sigma2; l];
    let mut lbf = vec![f64::NAN; l]; let mut kl = vec![f64::NAN; l];
    let mut q = vec![0.0; l * p];          // cached XtX b_l
    let mut xtxr = vec![0.0; p];
    let zm = (0..p).map(|j| { let bh = xty[j] / d[j]; let sh = (sigma2 / d[j]).sqrt(); (bh / sh).abs() }).filter(|x| !x.is_nan()).fold(0.0, f64::max);
    let mut elbo = vec![f64::NEG_INFINITY];
    let mut converged = false; let mut niter = 0;
    let mut xtr = vec![0.0; p]; let mut b = vec![0.0; p];
    let prof = std::env::var("FLASHLOCUS_PROFILE").is_ok();
    let (mut t_mv, mut t_ser, mut t_el) = (0.0f64, 0.0f64, 0.0f64);
    let t_start = std::time::Instant::now();
    for it in 1..=prm.max_iter {
        niter = it;
        for k in 0..l {
            let qk = &mut q[k * p..(k + 1) * p];
            for j in 0..p { xtxr[j] -= qk[j]; xtr[j] = xty[j] - xtxr[j] }
            let t0 = std::time::Instant::now();
            let s = ser(&xtr, &d, v[k], sigma2, &lpw, prm.estimate_prior_variance, prm.check_null_threshold);
            t_ser += t0.elapsed().as_secs_f64();
            let mut e1 = 0.0; let mut e2 = 0.0;
            for j in 0..p { let eb = s.alpha[j] * s.mu[j]; e1 += eb * xtr[j]; e2 += d[j] * s.alpha[j] * s.mu2[j] }
            kl[k] = -s.lbf_model + (-0.5 / sigma2 * (-2.0 * e1 + e2));
            for j in 0..p { b[j] = s.alpha[j] * s.mu[j] }
            let t0 = std::time::Instant::now();
            op.mul(&b, qk);
            t_mv += t0.elapsed().as_secs_f64();
            for j in 0..p { xtxr[j] += qk[j] }
            alpha[k * p..(k + 1) * p].copy_from_slice(&s.alpha);
            mu[k * p..(k + 1) * p].copy_from_slice(&s.mu);
            mu2[k * p..(k + 1) * p].copy_from_slice(&s.mu2);
            lbfv[k * p..(k + 1) * p].copy_from_slice(&s.lbf);
            v[k] = s.v; lbf[k] = s.lbf_model;
        }
        if prm.check_prior && v.iter().any(|x| *x > 100.0 * zm * zm) {
            return Err("estimated prior variance is unreasonably large (summary statistics / LD mismatch?)".into());
        }
        let t0 = std::time::Instant::now();
        // ELBO: Eloglik_ss - sum(KL), with ER2 from cached products
        let mut betabar = vec![0.0; p];
        let mut xb2 = 0.0; let mut postb2 = 0.0;
        for k in 0..l {
            for j in 0..p {
                let bk = alpha[k * p + j] * mu[k * p + j];
                betabar[j] += bk; xb2 += bk * q[k * p + j];
                postb2 += d[j] * alpha[k * p + j] * mu2[k * p + j];
            }
        }
        let er2 = yty - 2.0 * dot(&betabar, &xty) + dot(&betabar, &xtxr) - xb2 + postb2;
        let e = -n / 2.0 * (2.0 * std::f64::consts::PI * sigma2).ln() - 1.0 / (2.0 * sigma2) * er2 - kl.iter().sum::<f64>();
        if e.is_infinite() { return Err("objective became infinite".into()) }
        t_el += t0.elapsed().as_secs_f64();
        let prev = *elbo.last().unwrap();
        elbo.push(e);
        if e - prev < prm.tol { converged = true; break }
    }
    elbo.remove(0);
    let t0 = std::time::Instant::now();
    let sets = credible_sets(&alpha, &v, l, p, r, prm);
    if prof { eprintln!("profile p={p} iters={niter} total={:.4} matvec={t_mv:.4} ser={t_ser:.4} elbo={t_el:.4} cs={:.4}", t_start.elapsed().as_secs_f64() + t0.elapsed().as_secs_f64(), t0.elapsed().as_secs_f64()); }
    let mut pip = vec![0.0; p];
    let keep: Vec<usize> = (0..l).filter(|&k| v[k] > prm.prior_tol).collect();
    if !keep.is_empty() {
        for j in 0..p { let mut pr = 1.0; for &k in &keep { pr *= 1.0 - alpha[k * p + j] } pip[j] = 1.0 - pr }
    }
    Ok(Fit { p, l, alpha, mu, mu2, lbf_variable: lbfv, v, lbf, kl, elbo, niter, converged, pip, sets })
}

fn credible_sets<M: LdMat>(alpha: &[f64], v: &[f64], l: usize, p: usize, r: &M, prm: &Params) -> Vec<CredibleSet> {
    let mut cs_all: Vec<Vec<usize>> = Vec::with_capacity(l);
    let mut cov_all = vec![];
    for k in 0..l {
        let a = &alpha[k * p..(k + 1) * p];
        let mut o: Vec<usize> = (0..p).collect();
        o.sort_by(|i, j| a[*j].partial_cmp(&a[*i]).unwrap_or(std::cmp::Ordering::Equal)); // stable, decreasing
        let mut cum = 0.0; let mut nin = 0;
        for &j in &o { cum += a[j]; if cum < prm.coverage { nin += 1 } }
        let nin = nin + 1;
        let mut m: Vec<usize> = o[..nin.min(p)].to_vec(); m.sort();
        cov_all.push(m.iter().map(|&j| a[j]).sum::<f64>());
        cs_all.push(m);
    }
    let mut out = vec![];
    for k in 0..l {
        if !(v[k] > 1e-9) { continue }
        if cs_all[..k].iter().any(|c| *c == cs_all[k]) { continue }
        let m = &cs_all[k];
        let pur = if m.len() == 1 { [1.0, 1.0, 1.0] } else {
            // deterministic subsample of the first n_purity members by alpha when the set is large
            let pos: Vec<usize> = if m.len() > prm.n_purity {
                let a = &alpha[k * p..(k + 1) * p];
                let mut mm = m.clone(); mm.sort_by(|i, j| a[*j].partial_cmp(&a[*i]).unwrap()); mm.truncate(prm.n_purity); mm
            } else { m.clone() };
            let mut vals = vec![];
            for x in 0..pos.len() { for y in (x + 1)..pos.len() {
                let (i, j) = (pos[x], pos[y]);
                vals.push((r.get(i, j) / (r.get(i, i) * r.get(j, j)).sqrt()).abs());
            } }
            let mn = vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let h = vals.len();
            let med = if h % 2 == 1 { vals[h / 2] } else { (vals[h / 2 - 1] + vals[h / 2]) / 2.0 };
            [mn, mean, med]
        };
        if pur[0] >= prm.min_abs_corr { out.push(CredibleSet { effect: k, members: m.clone(), coverage: cov_all[k], purity: pur }) }
    }
    out.sort_by(|a, b| b.purity[0].partial_cmp(&a.purity[0]).unwrap());
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn susie_rss_matches_susier_0_14_2() {
        // tests/data/sim000: simulated locus (p 500, n 1e5) with the susieR 0.14.2 fit from the phase-1 benchmark
        let d = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data");
        let z: Vec<f64> = std::fs::read_to_string(format!("{d}/sim000.z")).unwrap().lines().map(|l| l.trim().parse().unwrap()).collect();
        let ld = crate::ldmat::read_full(&format!("{d}/sim000.ld.bin"), z.len()).unwrap();
        let f = super::susie_rss(&z, &ld, 100000.0, &super::Params::default()).unwrap();
        let pip: Vec<f64> = std::fs::read_to_string(format!("{d}/sim000.susier014.pip")).unwrap().lines().map(|l| l.trim().parse().unwrap()).collect();
        let m = f.pip.iter().zip(&pip).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(m < 1e-6, "max PIP difference {m}");
        let cs = std::fs::read_to_string(format!("{d}/sim000.susier014.cs")).unwrap();
        assert_eq!(f.sets.len(), cs.lines().filter(|l| !l.trim().is_empty()).count(), "number of credible sets");
    }
}
