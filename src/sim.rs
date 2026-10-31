//! Benchmark loci simulated from a PLINK reference (same design as PLAN.md phase 1):
//! contiguous block of p variants with MAF >= 1%, in-sample R, k in {1,2,3} causal variants,
//! h2 ~ U(0.0005, 0.002) split evenly, n = 100,000; z = sqrt(n) R[:,C] beta + X'g/sqrt(N).
use crate::{bfile, rng::Rng, simd};
use std::io::Write;

pub fn run(bfile: &str, out: &str, n_loci: usize, seed: u64, sizes: &[usize]) {
    std::fs::create_dir_all(out).unwrap();
    let b = bfile::Bfile::open(bfile).unwrap();
    let auto: Vec<usize> = (0..b.m()).filter(|&i| (1..=22).contains(&b.chr[i])).collect();
    let mut rng = Rng::new(seed);
    let mut truth = std::fs::File::create(format!("{out}/truth.tsv")).unwrap();
    writeln!(truth, "id\tp\tk\th2\tcausal\tchr\tstart\tn\tmax_abs_z").unwrap();
    let mut man = std::fs::File::create(format!("{out}/manifest.tsv")).unwrap();
    let mut manf = std::fs::File::create(format!("{out}/manifest_factor.tsv")).unwrap();
    let nn = b.n; let n_gwas = 100000.0f64;
    let mut l = 0;
    while l < n_loci {
        let p = sizes[l % sizes.len()];
        let s = auto[rng.below(auto.len() - p * 4)];
        if b.chr[s] != b.chr[s + p * 4 - 1] { continue }
        let mut g = b.dosage(s, s + p * 4);
        let keep: Vec<usize> = (0..p * 4).filter(|&i| {
            let r = &g[i * nn..(i + 1) * nn]; let af = r.iter().sum::<f64>() / (2.0 * nn as f64);
            let m = 2.0 * af; let var = r.iter().map(|x| (x - m) * (x - m)).sum::<f64>();
            af.min(1.0 - af) >= 0.01 && var > 0.0
        }).collect();
        if keep.len() < p { continue }
        let keep = &keep[..p];
        let mut x = Vec::with_capacity(p * nn);
        for &i in keep { x.extend_from_slice(&g[i * nn..(i + 1) * nn]) }
        g.clear();
        bfile::standardise(&mut x, nn);
        let r = bfile::ld(&x, p, nn);
        let k = 1 + l % 3;
        let h2 = 0.0005 + 0.0015 * rng.uniform();
        let mut c: Vec<usize> = vec![];
        while c.len() < k { let j = rng.below(p); if !c.contains(&j) { c.push(j) } }
        let beta: Vec<f64> = (0..k).map(|_| if rng.uniform() < 0.5 { -1.0 } else { 1.0 } * (h2 / k as f64).sqrt()).collect();
        let gn: Vec<f64> = (0..nn).map(|_| rng.normal()).collect();
        let z: Vec<f64> = (0..p).map(|j| {
            let sig: f64 = c.iter().zip(&beta).map(|(&cc, bb)| r[j * p + cc] * bb).sum::<f64>() * n_gwas.sqrt();
            sig + simd::dot(&x[j * nn..(j + 1) * nn], &gn) / (nn as f64).sqrt()
        }).collect();
        let id = format!("sim{l:03}");
        let mut zf = std::io::BufWriter::new(std::fs::File::create(format!("{out}/{id}.z")).unwrap());
        for v in &z { writeln!(zf, "{v:.8}").unwrap() }
        let bytes: Vec<u8> = r.iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(format!("{out}/{id}.ld.bin"), bytes).unwrap();
        std::fs::write(format!("{out}/{id}.geno.bin"), x.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
        let mut cs = c.clone(); cs.sort();
        writeln!(truth, "{id}\t{p}\t{k}\t{h2:.6}\t{}\t{}\t{}\t{}\t{:.3}", cs.iter().map(|x| (x + 1).to_string()).collect::<Vec<_>>().join(","),
                 b.chr[s], b.pos[s + keep[0]], n_gwas as u64, z.iter().fold(0.0f64, |a, v| a.max(v.abs()))).unwrap();
        writeln!(man, "{id}\t{out}/{id}.z\t{out}/{id}.ld.bin\t{}", n_gwas as u64).unwrap();
        writeln!(manf, "{id}\t{out}/{id}.z\t{out}/{id}.geno.bin:{nn}\t{}", n_gwas as u64).unwrap();
        eprintln!("{id} p={p} k={k}");
        l += 1;
    }
}
