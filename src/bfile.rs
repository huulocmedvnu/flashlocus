//! PLINK 1 binary reader (.bed SNP-major, .bim, .fam) and in-sample LD (correlation) with SIMD + threads.
use crate::simd;
use rayon::prelude::*;
use std::os::unix::fs::FileExt;

pub struct Bfile { pub n: usize, pub male: Vec<bool>, nb: usize, bed: std::fs::File, bim: Vec<u8>, line: Vec<u64>, pub chr: Vec<u8>, pub pos: Vec<u32> }

fn chr_code(c: &[u8]) -> u8 {
    let c = c.strip_prefix(b"chr").unwrap_or(c);
    match c { b"X" => 23, b"Y" => 24, b"XY" => 25, b"MT" | b"M" => 26, _ => std::str::from_utf8(c).ok().and_then(|s| s.parse().ok()).unwrap_or(0) }
}

impl Bfile {
    pub fn open(prefix: &str) -> Result<Bfile, String> {
        let fam = std::fs::read_to_string(format!("{prefix}.fam")).map_err(|e| format!("{prefix}.fam: {e}"))?;
        let male: Vec<bool> = fam.lines().filter(|l| !l.trim().is_empty()).map(|l| l.split_whitespace().nth(4) == Some("1")).collect();
        let n = male.len();
        let bim = std::fs::read(format!("{prefix}.bim")).map_err(|e| format!("{prefix}.bim: {e}"))?;
        // parsed in parallel over newline-aligned chunks
        let mut cuts = vec![0usize];
        let step = (bim.len() / 256).max(1 << 16);
        while *cuts.last().unwrap() < bim.len() {
            let mut b = (*cuts.last().unwrap() + step).min(bim.len());
            while b < bim.len() && bim[b - 1] != b'\n' { b += 1 }
            cuts.push(b);
        }
        let parts: Vec<(Vec<u64>, Vec<u8>, Vec<u32>)> = cuts.par_windows(2).map(|w| {
            let (mut line, mut chr, mut pos) = (vec![], vec![], vec![]);
            let mut st = w[0];
            while st < w[1] {
                let en = bim[st..w[1]].iter().position(|&c| c == b'\n').map(|x| st + x).unwrap_or(w[1]);
                let l = &bim[st..en];
                if !l.iter().all(|c| c.is_ascii_whitespace()) {
                    let mut f = l.split(|c| *c == b'\t' || *c == b' ').filter(|x| !x.is_empty());
                    let c = f.next().unwrap(); let _id = f.next(); let _cm = f.next(); let ps = f.next().unwrap();
                    line.push(st as u64); chr.push(chr_code(c));
                    pos.push(std::str::from_utf8(ps).unwrap().trim().parse().unwrap_or(0));
                }
                st = en + 1;
            }
            (line, chr, pos)
        }).collect();
        let (mut line, mut chr, mut pos) = (vec![], vec![], vec![]);
        for (l, c, p) in parts { line.extend(l); chr.extend(c); pos.extend(p) }
        let bed = std::fs::File::open(format!("{prefix}.bed")).map_err(|e| format!("{prefix}.bed: {e}"))?;
        let mut magic = [0u8; 3]; bed.read_exact_at(&mut magic, 0).map_err(|e| e.to_string())?;
        if magic != [0x6c, 0x1b, 0x01] { return Err("not a SNP-major PLINK .bed".into()) }
        Ok(Bfile { n, male, nb: (n + 3) / 4, bed, bim, line, chr, pos })
    }
    pub fn m(&self) -> usize { self.chr.len() }
    /// (id, a1, a2) of variant i
    pub fn info(&self, i: usize) -> (String, String, String) {
        let st = self.line[i] as usize;
        let en = self.bim[st..].iter().position(|c| *c == b'\n').map(|x| st + x).unwrap_or(self.bim.len());
        let f: Vec<&str> = std::str::from_utf8(&self.bim[st..en]).unwrap().split(|c| c == '\t' || c == ' ').filter(|x| !x.is_empty()).collect();
        (f[1].to_string(), f[4].to_string(), f[5].to_string())
    }
    /// A1 dosage (0/1/2), missing imputed by the variant mean; rows i0..i1. Returns row-major (i1-i0) x n.
    pub fn dosage(&self, i0: usize, i1: usize) -> Vec<f64> {
        const LUT: [f64; 4] = [2.0, f64::NAN, 1.0, 0.0];
        let mut raw = vec![0u8; (i1 - i0) * self.nb];
        self.bed.read_exact_at(&mut raw, 3 + (i0 * self.nb) as u64).expect("bed read");
        let n = self.n;
        let mut g = vec![0.0f64; (i1 - i0) * n];
        g.par_chunks_mut(n).zip(raw.par_chunks(self.nb)).for_each(|(row, bytes)| {
            let (mut s, mut c) = (0.0, 0usize);
            for k in 0..n { let v = LUT[((bytes[k / 4] >> (2 * (k % 4))) & 3) as usize]; row[k] = v; if !v.is_nan() { s += v; c += 1 } }
            let m = if c > 0 { s / c as f64 } else { 0.0 };
            for v in row.iter_mut() { if v.is_nan() { *v = m } }
        });
        g
    }
}

/// Standardise rows (mean 0, population sd 1) in place; returns allele frequency of A1 and sd per row.
pub fn standardise(g: &mut [f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let stats: Vec<(f64, f64)> = g.par_chunks_mut(n).map(|r| {
        let m = r.iter().sum::<f64>() / n as f64;
        let sd = (r.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64).sqrt();
        for x in r.iter_mut() { *x = if sd > 0.0 { (*x - m) / sd } else { 0.0 } }
        (m / 2.0, sd)
    }).collect();
    (stats.iter().map(|s| s.0).collect(), stats.iter().map(|s| s.1).collect())
}

/// R = X X' / n for standardised rows X (p x n). Upper triangle by rows in parallel, then mirrored.
pub fn ld(x: &[f64], p: usize, n: usize) -> Vec<f64> {
    let mut r = vec![0.0f64; p * p];
    r.par_chunks_mut(p).enumerate().for_each(|(i, row)| {
        let xi = &x[i * n..(i + 1) * n];
        for j in i..p { row[j] = simd::dot(xi, &x[j * n..(j + 1) * n]) / n as f64 }
        row[i] = 1.0;
    });
    for i in 0..p { for j in 0..i { r[i * p + j] = r[j * p + i] } }
    r
}

/// Exact LD from 2-bit genotypes when the window has no missing calls: dosage x = a + b with bit planes
/// a = [x >= 1], b = [x == 2], so sum xy = |a1&a2| + |a1&b2| + |b1&a2| + |b1&b2| and sum x^2 = |a| + 3|b|.
/// r = (n sxy - sx sy) / sqrt((n sxx - sx^2)(n syy - sy^2)), integer sums as in PLINK. Returns None if any call
/// is missing (the caller falls back to mean-imputed floating point). Output: kept indices, A1 freq, full p x p.
impl Bfile {
    /// x_haploid: on chrX (code 23) male calls are haploid as in PLINK 1.9 (hom A1 = 1, hom A2 = 0, het = missing).
    pub fn ld_exact(&self, idx: &[usize], maf: f64, x_haploid: bool) -> Option<(Vec<usize>, Vec<f64>, Vec<f64>)> {
        let n = self.n; let w = (n + 63) / 64; let nb = self.nb;
        let i0 = idx[0]; let i1 = idx[idx.len() - 1] + 1;
        let mut raw = vec![0u8; (i1 - i0) * nb];
        self.bed.read_exact_at(&mut raw, 3 + (i0 * nb) as u64).expect("bed read");
        // planes per variant: [a words..., b words...]
        let planes: Vec<Option<(Vec<u64>, i64, i64)>> = idx.par_iter().map(|&i| {
            let by = &raw[(i - i0) * nb..(i - i0 + 1) * nb];
            let mut pl = vec![0u64; 2 * w];
            let hap = x_haploid && self.chr[i] == 23;
            for k in 0..n {
                let c = (by[k / 4] >> (2 * (k % 4))) & 3;
                let hm = hap && self.male[k];
                match c {
                    0 => { pl[k / 64] |= 1 << (k % 64); if !hm { pl[w + k / 64] |= 1 << (k % 64) } }
                    2 => { if hm { return None } pl[k / 64] |= 1 << (k % 64) }
                    3 => {}
                    _ => return None,
                }
            }
            let pa: i64 = pl[..w].iter().map(|x| x.count_ones() as i64).sum();
            let pb: i64 = pl[w..].iter().map(|x| x.count_ones() as i64).sum();
            Some((pl, pa + pb, pa + 3 * pb))
        }).collect();
        if planes.iter().any(|x| x.is_none()) { return None }
        let planes: Vec<(Vec<u64>, i64, i64)> = planes.into_iter().map(|x| x.unwrap()).collect();
        let ni = n as i64;
        let keep: Vec<usize> = (0..idx.len()).filter(|&k| {
            let (sx, sxx) = (planes[k].1, planes[k].2);
            let af = sx as f64 / self.allele_count(idx[k], x_haploid);
            af.min(1.0 - af) >= maf && ni * sxx - sx * sx > 0
        }).collect();
        let p = keep.len();
        let flat: Vec<u64> = keep.iter().flat_map(|&k| planes[k].0.iter().cloned()).collect();
        let sx: Vec<i64> = keep.iter().map(|&k| planes[k].1).collect();
        let dx: Vec<f64> = keep.iter().map(|&k| ((ni * planes[k].2 - planes[k].1 * planes[k].1) as f64).sqrt()).collect();
        let mut r = vec![0.0f64; p * p];
        r.par_chunks_mut(p).enumerate().for_each(|(i, row)| {
            let (ai, bi) = (&flat[i * 2 * w..i * 2 * w + w], &flat[i * 2 * w + w..(i + 1) * 2 * w]);
            for j in 0..p {
                let (aj, bj) = (&flat[j * 2 * w..j * 2 * w + w], &flat[j * 2 * w + w..(j + 1) * 2 * w]);
                let mut sxy = 0i64;
                for t in 0..w { sxy += ((ai[t] & aj[t]).count_ones() + (ai[t] & bj[t]).count_ones() + (bi[t] & aj[t]).count_ones() + (bi[t] & bj[t]).count_ones()) as i64 }
                row[j] = (ni * sxy - sx[i] * sx[j]) as f64 / (dx[i] * dx[j]);
            }
            row[i] = 1.0;
        });
        let af = keep.iter().map(|&k| planes[k].1 as f64 / self.allele_count(idx[k], x_haploid)).collect();
        Some((keep.iter().map(|&k| idx[k]).collect(), af, r))
    }
}

impl Bfile {
    /// Number of allele copies at variant i (2n, or 2*females + males on chrX with haploid males).
    pub fn allele_count(&self, i: usize, x_haploid: bool) -> f64 {
        if x_haploid && self.chr[i] == 23 { let m = self.male.iter().filter(|&&x| x).count(); (2 * (self.n - m) + m) as f64 } else { 2.0 * self.n as f64 }
    }
    /// Float-path dosage with PLINK 1.9 chrX coding: male calls halved, male hets missing (then mean-imputed).
    pub fn dosage_x(&self, idx: &[usize], x_haploid: bool) -> Vec<f64> {
        let i0 = idx[0];
        let gr = self.dosage(i0, idx[idx.len() - 1] + 1);
        let n = self.n;
        let mut g = Vec::with_capacity(idx.len() * n);
        for &i in idx {
            let row = &gr[(i - i0) * n..(i - i0 + 1) * n];
            if x_haploid && self.chr[i] == 23 {
                let mut r: Vec<f64> = row.iter().zip(&self.male).map(|(&v, &m)| if m { if v == 1.0 { f64::NAN } else { v / 2.0 } } else { v }).collect();
                let (s, c) = r.iter().filter(|v| !v.is_nan()).fold((0.0, 0usize), |a, v| (a.0 + v, a.1 + 1));
                let mu = if c > 0 { s / c as f64 } else { 0.0 };
                for v in r.iter_mut() { if v.is_nan() { *v = mu } }
                g.extend(r);
            } else { g.extend_from_slice(row) }
        }
        g
    }
}

impl Bfile {
    /// LD for variants on `chr` in [start, end] (optionally only IDs in `extract`), MAF >= maf.
    /// Exact integer path unless a call is missing or `force_float`. Returns (bim indices kept, A1 freq, p x p r, mode).
    pub fn ld_region(&self, chr: u8, start: u32, end: u32, extract: Option<&std::collections::HashSet<String>>, maf: f64,
                     x_haploid: bool, force_float: bool) -> Option<(Vec<usize>, Vec<f64>, Vec<f64>, &'static str)> {
        let all: Vec<usize> = (0..self.m()).filter(|&i| self.chr[i] == chr && self.pos[i] >= start && self.pos[i] <= end).collect();
        let idx: Vec<usize> = match extract { Some(x) => all.into_iter().filter(|&i| x.contains(&self.info(i).0)).collect(), None => all };
        if idx.is_empty() { return None }
        if !force_float { if let Some((k, af, r)) = self.ld_exact(&idx, maf, x_haploid) { return Some((k, af, r, "exact-int")) } }
        let mut g = self.dosage_x(&idx, x_haploid);
        let (mut af, sd) = standardise(&mut g, self.n);
        // standardise() returns mean/2; convert to A1 frequency over allele copies on haploid-male chrX
        for (k, &i) in idx.iter().enumerate() { if x_haploid && self.chr[i] == 23 { af[k] = af[k] * 2.0 * self.n as f64 / self.allele_count(i, true) } }
        let kk: Vec<usize> = (0..idx.len()).filter(|&k| af[k].min(1.0 - af[k]) >= maf && sd[k] > 0.0).collect();
        let mut x = Vec::with_capacity(kk.len() * self.n);
        for &k in &kk { x.extend_from_slice(&g[k * self.n..(k + 1) * self.n]) }
        let r = ld(&x, kk.len(), self.n);
        Some((kk.iter().map(|&k| idx[k]).collect(), kk.iter().map(|&k| af[k]).collect(), r, "float"))
    }
}

/// Write f64 slice as little-endian bytes without an intermediate copy.
pub fn write_f64(path: &str, v: &[f64]) {
    #[cfg(target_endian = "little")]
    {
        // SAFETY: f64 has no padding; reinterpreting the slice as bytes is valid on little-endian targets.
        let b = unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, v.len() * 8) };
        std::fs::write(path, b).unwrap_or_else(|e| panic!("{path}: {e}"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tiny PLINK fileset: n samples (alternating sex), m variants on `chr`, genotype codes from a fixed LCG.
    fn write_bfile(dir: &std::path::Path, chr: u8, n: usize, m: usize, missing: bool) -> String {
        let pre = dir.join(format!("t{chr}{missing}")).to_string_lossy().to_string();
        let fam: String = (0..n).map(|k| format!("f{k} i{k} 0 0 {} -9\n", 1 + k % 2)).collect();
        std::fs::write(format!("{pre}.fam"), fam).unwrap();
        let bim: String = (0..m).map(|j| format!("{chr}\trs{j}\t0\t{}\tA\tG\n", 1000 + j)).collect();
        std::fs::write(format!("{pre}.bim"), bim).unwrap();
        let nb = (n + 3) / 4;
        let mut bed = vec![0x6c, 0x1b, 0x01];
        let mut s: u64 = 20261031;
        for j in 0..m {
            let mut row = vec![0u8; nb];
            for k in 0..n {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let male_x = chr == 23 && k % 2 == 0;
                // codes: 0 hom A1, 2 het, 3 hom A2, 1 missing; male chrX calls are homozygous
                let mut c = match (s >> 33) % 4 { 0 => 0u8, 1 | 2 => if male_x { 3 } else { 2 }, _ => 3 };
                if missing && j == 3 && k == 5 { c = 1 }
                row[k / 4] |= c << (2 * (k % 4));
            }
            bed.extend(row);
        }
        std::fs::write(format!("{pre}.bed"), bed).unwrap();
        pre
    }

    fn float_ld(b: &Bfile, idx: &[usize], xh: bool) -> Vec<f64> {
        let mut g = b.dosage_x(idx, xh);
        standardise(&mut g, b.n);
        ld(&g, idx.len(), b.n)
    }

    #[test]
    fn exact_ld_equals_float_ld_autosome_and_chrx() {
        let dir = std::env::temp_dir().join(format!("flashlocus_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for chr in [1u8, 23] {
            let b = Bfile::open(&write_bfile(&dir, chr, 203, 40, false)).unwrap();
            let idx: Vec<usize> = (0..b.m()).collect();
            for xh in [true, false] {
                let (keep, _, r) = b.ld_exact(&idx, 0.0, xh).expect("no missing calls");
                assert_eq!(keep.len(), idx.len());
                let f = float_ld(&b, &idx, xh);
                let m = r.iter().zip(&f).map(|(a, c)| (a - c).abs()).fold(0.0, f64::max);
                assert!(m < 1e-12, "chr {chr} x_haploid {xh}: max diff {m}");
            }
        }
        // a missing call switches to the floating-point path
        let b = Bfile::open(&write_bfile(&dir, 1, 203, 40, true)).unwrap();
        assert!(b.ld_exact(&(0..b.m()).collect::<Vec<_>>(), 0.0, true).is_none());
        std::fs::remove_dir_all(&dir).ok();
    }
}
