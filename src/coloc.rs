//! Exact port of coloc 5.2.3 coloc.bf_bf (as called by coloc.susie on credible-set rows of lbf_variable).
//! Kept conventions of the R code, including the ones that look odd:
//!  - logbf_to_pp(last_is_null = TRUE) treats the LAST column as the null when no "null" column exists,
//!    so the last SNP gets prior 1 - n*pi and every row is shifted by its last-column value;
//!  - trimming uses the full matrices, combine.abf uses the raw lbf restricted to shared SNPs;
//!  - shared SNPs follow the column order of bf1 (R intersect), hits are the first maximum.

pub struct Pair { pub i: usize, pub j: usize, pub nsnps: usize, pub hit1: usize, pub hit2: usize, pub pp: [f64; 5] }

#[allow(dead_code)] // nsnps/all_dropped are part of the API, not used by the atlas writer
pub struct Out { pub nsnps: usize, pub pairs: Vec<Pair>, pub all_dropped: bool, pub bf2_unaligned: bool }

#[inline]
fn logsum(x: &[f64]) -> f64 {
    let m = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + x.iter().map(|v| (v - m).exp()).sum::<f64>().ln()
}

fn logdiff(x: f64, y: f64) -> f64 { let m = x.max(y); m + ((x - m).exp() - (y - m).exp()).ln() }

/// pp row of logbf_to_pp(bf, pi, last_is_null = TRUE) for one lbf row.
fn pp_row(bf: &[f64], pi0: f64) -> Vec<f64> {
    let p = bf.len();
    let n = p - 1;
    let pi = if pi0 > 1.0 / n as f64 { 1.0 / n as f64 } else { pi0 };
    let mut pri: Vec<f64> = (0..p).map(|k| if k < n { pi } else { 1.0 - n as f64 * pi }).collect();
    if pri.iter().any(|&v| v == 0.0) {
        for v in pri.iter_mut() { if *v == 0.0 { *v = 1e-16 } }
        let s: f64 = pri.iter().sum();
        for v in pri.iter_mut() { *v /= s }
    }
    let last = bf[p - 1];
    let t: Vec<f64> = bf.iter().zip(&pri).map(|(b, q)| (b - last) + q.ln()).collect();
    let d = logsum(&t);
    t.iter().map(|v| (v - d).exp()).collect()
}

pub fn combine_abf(l1: &[f64], l2: &[f64], p1: f64, p2: f64, p12: f64) -> [f64; 5] {
    let lsum: Vec<f64> = l1.iter().zip(l2).map(|(a, b)| a + b).collect();
    let (s1, s2, s12) = (logsum(l1), logsum(l2), logsum(&lsum));
    let a = [0.0, p1.ln() + s1, p2.ln() + s2, p1.ln() + p2.ln() + logdiff(s1 + s2, s12), p12.ln() + s12];
    let d = logsum(&a);
    [(a[0] - d).exp(), (a[1] - d).exp(), (a[2] - d).exp(), (a[3] - d).exp(), (a[4] - d).exp()]
}

fn which_max(x: &[f64]) -> usize {
    let mut b = 0;
    for k in 1..x.len() { if x[k] > x[b] { b = k } }  // R which.max: first maximum, NA skipped
    b
}

/// bf1: rows of lbf over keys1 (p1 columns); bf2 likewise. Returns pairs (row i of bf1, row j of bf2).
/// hit1/hit2 are column indices into keys1/keys2 (the SNP names R reports).
pub fn bf_bf<K: std::hash::Hash + Eq>(bf1: &[&[f64]], keys1: &[K], bf2: &[&[f64]], keys2: &[K],
                                       p1: f64, p2: f64, p12: f64, overlap_min: f64) -> Out {
    let pos2: std::collections::HashMap<&K, usize> = {
        let mut m = std::collections::HashMap::with_capacity(keys2.len());
        for (k, key) in keys2.iter().enumerate() { m.entry(key).or_insert(k); }  // match(): first occurrence
        m
    };
    let mut seen = std::collections::HashSet::new();
    let mut c1 = Vec::new(); let mut c2 = Vec::new();  // shared SNPs: column in bf1, column in bf2
    for (k, key) in keys1.iter().enumerate() {
        if let Some(&j) = pos2.get(key) { if seen.insert(key) { c1.push(k); c2.push(j) } }
    }
    let ns = c1.len();
    if ns == 0 { return Out { nsnps: 0, pairs: vec![], all_dropped: false, bf2_unaligned: false } }
    let prop = |bf: &[f64], cols: &[usize], p: f64| -> f64 {
        let pp = pp_row(bf, p);
        cols.iter().map(|&k| pp[k]).sum::<f64>() / pp.iter().sum::<f64>()
    };
    let prop1: Vec<f64> = bf1.iter().map(|r| prop(r, &c1, p1)).collect();
    let prop2: Vec<f64> = bf2.iter().map(|r| prop(r, &c2, p2)).collect();
    // expand.grid(i, j): i varies fastest
    let mut todo = Vec::new();
    for j in 0..bf2.len() { for i in 0..bf1.len() { todo.push((i, j)) } }
    let keep: Vec<(usize, usize)> = todo.iter().cloned().filter(|&(i, j)| !(prop1[i] < overlap_min || prop2[j] < overlap_min)).collect();
    if keep.is_empty() {
        let pairs = todo.iter().map(|&(i, j)| {
            let h0 = (1.0 - pp_row(bf1[i], p1).iter().sum::<f64>()).min(1.0 - pp_row(bf2[j], p2).iter().sum::<f64>());
            Pair { i, j, nsnps: ns, hit1: which_max(&pp_row(bf1[i], p1)), hit2: which_max(&pp_row(bf2[j], p2)), pp: [h0, f64::NAN, f64::NAN, f64::NAN, f64::NAN] }
        }).collect();
        return Out { nsnps: ns, pairs, all_dropped: true, bf2_unaligned: false };
    }
    // R subsets bf2 to match(isnps, colnames(bf2)) only when some bf2 SNPs are not shared; otherwise bf2
    // keeps its own column order and is paired with bf1 by position. Reproduced; flagged when orders differ.
    let p2n = bf2[0].len();
    let c2u: Vec<usize> = if ns == p2n { (0..p2n).collect() } else { c2.clone() };
    let bf2_unaligned = c2u != c2;
    let pairs = keep.iter().map(|&(i, j)| {
        let l1: Vec<f64> = c1.iter().map(|&k| bf1[i][k]).collect();
        let l2: Vec<f64> = c2u.iter().map(|&k| bf2[j][k]).collect();
        let mut pp = combine_abf(&l1, &l2, p1, p2, p12);
        let s: Vec<f64> = l1.iter().zip(&l2).map(|(a, b)| a + b).collect();
        let d = logsum(&s);
        if s.iter().all(|v| (v - d).exp().is_nan()) { pp = [1.0, 0.0, 0.0, 0.0, 0.0] }
        Pair { i, j, nsnps: ns, hit1: c1[which_max(&l1)], hit2: c2u[which_max(&l2)], pp }
    }).collect();
    Out { nsnps: ns, pairs, all_dropped: false, bf2_unaligned }
}

#[cfg(test)]
mod tests {
    // Fixtures written by coloc 5.2.3 coloc.bf_bf (tests/data/coloc_*): bf matrices stored SNP x row.
    fn read_bf(p: &str) -> (Vec<String>, Vec<Vec<f64>>) {
        let t = std::fs::read_to_string(p).unwrap();
        let mut keys = vec![]; let mut cols: Vec<Vec<f64>> = vec![];
        for l in t.lines().skip(1) {
            let f: Vec<&str> = l.split('\t').collect();
            keys.push(f[0].to_string());
            if cols.is_empty() { cols = vec![vec![]; f.len() - 1] }
            for (k, v) in f[1..].iter().enumerate() { cols[k].push(v.parse().unwrap()) }
        }
        (keys, cols)
    }

    fn check(case: &str) {
        let d = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data");
        let (k1, b1) = read_bf(&format!("{d}/coloc_{case}_bf1.tsv"));
        let (k2, b2) = read_bf(&format!("{d}/coloc_{case}_bf2.tsv"));
        let r1: Vec<&[f64]> = b1.iter().map(|r| r.as_slice()).collect();
        let r2: Vec<&[f64]> = b2.iter().map(|r| r.as_slice()).collect();
        let o = super::bf_bf(&r1, &k1, &r2, &k2, 1e-4, 1e-4, 5e-6, 0.5);
        let exp = std::fs::read_to_string(format!("{d}/coloc_{case}_R.tsv")).unwrap();
        let rows: Vec<Vec<String>> = exp.lines().skip(1).map(|l| l.split('\t').map(|s| s.trim_matches('"').to_string()).collect()).collect();
        assert_eq!(o.pairs.len(), rows.len(), "{case}: number of pairs");
        for r in &rows {
            let (i, j): (usize, usize) = (r[0].parse::<usize>().unwrap() - 1, r[1].parse::<usize>().unwrap() - 1);
            let p = o.pairs.iter().find(|p| p.i == i && p.j == j).expect("pair present");
            assert_eq!(p.nsnps.to_string(), r[2]);
            assert_eq!(k1[p.hit1], r[3]); assert_eq!(k2[p.hit2], r[4]);
            for h in 0..5 {
                let e: f64 = r[5 + h].parse().unwrap();
                if e.is_nan() { assert!(p.pp[h].is_nan(), "{case} H{h}") } else { assert!((p.pp[h] - e).abs() < 1e-10, "{case} H{h}: {} vs {e}", p.pp[h]) }
            }
        }
    }

    #[test] fn coloc_partial_overlap_matches_r() { check("partial") }
    #[test] fn coloc_all_trimmed_matches_r() { check("dropall") }
}
