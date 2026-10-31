//! Parity and speed against a reference tool (PLAN.md phase 1 gate).
//! Rust results: <id>.pip, <id>.cs.tsv; reference: <id>.pip, <id>.cs (effect<TAB>members, 1-based).
use std::collections::{BTreeSet, HashMap};

fn sets(path: &str, rust: bool) -> BTreeSet<Vec<usize>> {
    let Ok(t) = std::fs::read_to_string(path) else { return BTreeSet::new() };
    t.lines().skip(if rust { 1 } else { 0 }).filter(|l| !l.trim().is_empty()).map(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        let m = if rust { f[6] } else { f[1] };
        let mut v: Vec<usize> = m.split(',').filter_map(|x| x.trim().parse().ok()).collect(); v.sort(); v
    }).collect()
}
fn read_vec(p: &str) -> Option<Vec<f64>> { std::fs::read_to_string(p).ok().map(|t| t.lines().filter_map(|l| l.trim().parse().ok()).collect()) }
fn timing(p: &str) -> HashMap<String, f64> {
    let t = std::fs::read_to_string(p).unwrap_or_default(); let mut it = t.lines();
    let h: Vec<&str> = it.next().unwrap_or("").split('\t').collect();
    let (ci, cs) = (h.iter().position(|x| *x == "id").unwrap_or(0), h.iter().position(|x| *x == "seconds").unwrap_or(2));
    it.filter_map(|l| { let f: Vec<&str> = l.split('\t').collect(); f.get(cs)?.parse().ok().map(|s| (f[ci].to_string(), s)) }).collect()
}
fn median(mut v: Vec<f64>) -> f64 { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); let n = v.len(); if n == 0 { f64::NAN } else if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 } }

pub fn run(manifest: &str, rust_dir: &str, ref_dir: &str) {
    let ids: Vec<String> = std::fs::read_to_string(manifest).unwrap().lines().filter(|l| !l.is_empty()).map(|l| l.split('\t').next().unwrap().to_string()).collect();
    let (tr, tf) = (timing(&format!("{rust_dir}/timing.tsv")), timing(&format!("{ref_dir}/timing.tsv")));
    let (mut a, mut b) = (vec![], vec![]); let (mut same, mut n, mut maxd, mut big) = (0, 0, 0.0f64, 0); let mut ratio = vec![];
    for id in &ids {
        let (Some(x), Some(y)) = (read_vec(&format!("{rust_dir}/{id}.pip")), read_vec(&format!("{ref_dir}/{id}.pip"))) else { continue };
        if x.len() != y.len() { continue }
        let d = x.iter().zip(&y).map(|(u, v)| (u - v).abs()).fold(0.0, f64::max);
        maxd = maxd.max(d); if d > 0.01 { big += 1 }
        a.extend(x); b.extend(y); n += 1;
        if sets(&format!("{rust_dir}/{id}.cs.tsv"), true) == sets(&format!("{ref_dir}/{id}.cs"), false) { same += 1 }
        if let (Some(s1), Some(s2)) = (tr.get(id), tf.get(id)) { if *s1 > 0.0 { ratio.push(s2 / s1) } }
    }
    let m = a.len() as f64; let (ma, mb) = (a.iter().sum::<f64>() / m, b.iter().sum::<f64>() / m);
    let cov: f64 = a.iter().zip(&b).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let r = cov / (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>().sqrt() * b.iter().map(|y| (y - mb).powi(2)).sum::<f64>().sqrt());
    let tot = |t: &HashMap<String, f64>| ids.iter().filter_map(|i| t.get(i)).sum::<f64>();
    println!("reference\t{ref_dir}\nloci\t{n}\npip_pearson_r\t{r:.6}\nmax_abs_pip_diff\t{maxd:.4e}\nloci_pip_diff_gt_0.01\t{big}\ncs_identical_share\t{:.4}\nmedian_time_ratio_ref_over_rust\t{:.2}\ntotal_time_ratio\t{:.2}",
             same as f64 / n.max(1) as f64, median(ratio), tot(&tf) / tot(&tr));
}
