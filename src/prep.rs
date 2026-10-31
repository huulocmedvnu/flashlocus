//! Real-data loci: GWAS summary-statistic windows -> allele-aligned z, standardised reference genotypes
//! (factor), full LD (for reference tools). The reference .bim is scanned once for all windows.
//! windows list: id<TAB>path ; --format ssf | finngen | anm ; --match rsid | chrpos ; --n N (or per-SNP N column)
use crate::bfile::{self, Bfile};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::os::unix::fs::FileExt;

struct Row { key: String, ea: String, oa: String, z: f64, n: f64 }

fn cols(fmt: &str) -> (&'static str, &'static str, &'static str, &'static str, &'static str, &'static str, &'static str, &'static str) {
    // (rsid, chr, pos, effect allele, other allele, beta, se, n)
    match fmt {
        "finngen" => ("rsids", "#chrom", "pos", "alt", "ref", "beta", "sebeta", ""),
        "anm" => ("", "CHR", "POS", "Effect_Allele", "Other_Allele", "Effect", "SE", "N"),
        _ => ("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", ""),
    }
}

fn read_window(path: &str, fmt: &str, by_rsid: bool) -> Vec<Row> {
    let t = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut it = t.lines();
    let h: Vec<&str> = it.next().unwrap().split('\t').collect();
    let (crs, cch, cps, cea, coa, cb, cse, cn) = cols(fmt);
    let ix = |c: &str| -> Option<usize> { if c.is_empty() { None } else { h.iter().position(|x| *x == c).or_else(|| if c == "rsid" { h.iter().position(|x| *x == "rs_id") } else { None }) } };
    let (irs, ich, ips, iea, ioa, ib, ise, inn) = (ix(crs), ix(cch).unwrap(), ix(cps).unwrap(), ix(cea).unwrap(), ix(coa).unwrap(), ix(cb).unwrap(), ix(cse).unwrap(), ix(cn));
    it.filter_map(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        let (b, s): (f64, f64) = (f[ib].parse().ok()?, f[ise].parse().ok()?);
        if !(s > 0.0) || !b.is_finite() { return None }
        let key = if by_rsid { f[irs?].split(',').next()?.to_string() } else { format!("{}:{}", f[ich].trim_start_matches("chr"), f[ips]) };
        if by_rsid && !key.starts_with("rs") { return None }
        Some(Row { key, ea: f[iea].to_uppercase(), oa: f[ioa].to_uppercase(), z: b / s, n: inn.and_then(|i| f[i].parse().ok()).unwrap_or(f64::NAN) })
    }).collect()
}

fn comp(a: &str) -> String { a.chars().rev().map(|c| match c { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', x => x }).collect() }

pub fn run(bfile: &str, list: &str, fmt: &str, by_rsid: bool, n_fixed: Option<f64>, maf: f64, out: &str) {
    std::fs::create_dir_all(out).unwrap();
    let b = Bfile::open(bfile).unwrap();
    let wins: Vec<(String, String)> = std::fs::read_to_string(list).unwrap().lines().filter(|l| !l.is_empty()).map(|l| { let f: Vec<&str> = l.split('\t').collect(); (f[0].to_string(), f[1].to_string()) }).collect();
    let rows: Vec<Vec<Row>> = wins.par_iter().map(|(_, p)| read_window(p, fmt, by_rsid)).collect();
    // key -> bim indices (all needed keys, one pass over the reference)
    let need: std::collections::HashSet<&str> = rows.iter().flatten().map(|r| r.key.as_str()).collect();
    let keys: Vec<Option<String>> = (0..b.m()).into_par_iter().map(|i| {
        let k = if by_rsid { b.info(i).0 } else { format!("{}:{}", b.chr[i], b.pos[i]) };
        if need.contains(k.as_str()) { Some(k) } else { None }
    }).collect();
    let mut map: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, k) in keys.into_iter().enumerate() { if let Some(k) = k { map.entry(k).or_default().push(i) } }
    eprintln!("prep: {} windows, {} keys needed, {} found in reference", wins.len(), need.len(), map.len());
    let nb = (b.n + 3) / 4;
    let bed = std::fs::File::open(format!("{bfile}.bed")).unwrap();
    let mut man = std::fs::File::create(format!("{out}/manifest.tsv")).unwrap();
    let mut manf = std::fs::File::create(format!("{out}/manifest_factor.tsv")).unwrap();
    let mut summ = std::fs::File::create(format!("{out}/prep_summary.tsv")).unwrap();
    writeln!(summ, "id\tin_window\tmatched\tflipped\tdropped_palindromic\tdropped_allele\tdropped_maf\tp\tn\tmax_abs_z").unwrap();
    let results: Vec<String> = wins.par_iter().zip(rows.par_iter()).map(|((id, _), rs)| {
        let (mut flip, mut pal, mut bad) = (0, 0, 0);
        let mut sel: Vec<(usize, f64, f64)> = vec![];
        for r in rs {
            let Some(cands) = map.get(&r.key) else { continue };
            let mut hit = None;
            for &i in cands {
                let (_, a1, a2) = b.info(i); let (a1, a2) = (a1.to_uppercase(), a2.to_uppercase());
                let palin = comp(&a1) == a2;
                if palin { pal += 1; hit = Some(None); break }
                if r.ea == a1 && r.oa == a2 || comp(&r.ea) == a1 && comp(&r.oa) == a2 { hit = Some(Some((i, r.z))); break }
                if r.ea == a2 && r.oa == a1 || comp(&r.ea) == a2 && comp(&r.oa) == a1 { flip += 1; hit = Some(Some((i, -r.z))); break }
            }
            match hit { Some(Some((i, z))) => sel.push((i, z, r.n)), Some(None) => {}, None => bad += 1 }
        }
        sel.sort_by_key(|x| x.0); sel.dedup_by_key(|x| x.0);
        let mut g = vec![0.0f64; sel.len() * b.n]; let mut raw = vec![0u8; nb];
        const LUT: [f64; 4] = [2.0, f64::NAN, 1.0, 0.0];
        for (k, (i, _, _)) in sel.iter().enumerate() {
            bed.read_exact_at(&mut raw, 3 + (*i * nb) as u64).unwrap();
            let row = &mut g[k * b.n..(k + 1) * b.n]; let (mut s, mut c) = (0.0, 0);
            for j in 0..b.n { let v = LUT[((raw[j / 4] >> (2 * (j % 4))) & 3) as usize]; row[j] = v; if !v.is_nan() { s += v; c += 1 } }
            let m = if c > 0 { s / c as f64 } else { 0.0 }; for v in row.iter_mut() { if v.is_nan() { *v = m } }
        }
        let (af, sd) = bfile::standardise(&mut g, b.n);
        let keep: Vec<usize> = (0..sel.len()).filter(|&k| af[k].min(1.0 - af[k]) >= maf && sd[k] > 0.0).collect();
        let dmaf = sel.len() - keep.len();
        let p = keep.len();
        let mut x = Vec::with_capacity(p * b.n); for &k in &keep { x.extend_from_slice(&g[k * b.n..(k + 1) * b.n]) }
        let zs: Vec<f64> = keep.iter().map(|&k| sel[k].1).collect();
        let mut ns: Vec<f64> = keep.iter().map(|&k| sel[k].2).filter(|v| v.is_finite()).collect();
        ns.sort_by(|a, c| a.partial_cmp(c).unwrap());
        let n = n_fixed.unwrap_or_else(|| if ns.is_empty() { f64::NAN } else { ns[ns.len() / 2] });
        let mut zf = std::io::BufWriter::new(std::fs::File::create(format!("{out}/{id}.z")).unwrap());
        for v in &zs { writeln!(zf, "{v:.8}").unwrap() }
        std::fs::write(format!("{out}/{id}.geno.bin"), x.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
        let r = bfile::ld(&x, p, b.n);
        std::fs::write(format!("{out}/{id}.ld.bin"), r.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>()).unwrap();
        let mut sf = std::io::BufWriter::new(std::fs::File::create(format!("{out}/{id}.snps.tsv")).unwrap());
        writeln!(sf, "id\tchr\tpos\ta1\ta2\taf_a1\tz").unwrap();
        for (q, &k) in keep.iter().enumerate() { let i = sel[k].0; let (sid, a1, a2) = b.info(i); writeln!(sf, "{sid}\t{}\t{}\t{a1}\t{a2}\t{:.5}\t{:.5}", b.chr[i], b.pos[i], af[k], zs[q]).unwrap() }
        let mz = zs.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        format!("{id}\t{}\t{}\t{flip}\t{pal}\t{bad}\t{dmaf}\t{p}\t{n}\t{mz:.2}", rs.len(), sel.len())
    }).collect();
    for (k, l) in results.iter().enumerate() {
        writeln!(summ, "{l}").unwrap();
        let f: Vec<&str> = l.split('\t').collect(); let id = &wins[k].0;
        if f[7].parse::<usize>().unwrap_or(0) >= 10 && f[8] != "NaN" {
            writeln!(man, "{id}\t{out}/{id}.z\t{out}/{id}.ld.bin\t{}", f[8]).unwrap();
            writeln!(manf, "{id}\t{out}/{id}.z\t{out}/{id}.geno.bin:{}\t{}", b.n, f[8]).unwrap();
        }
    }
}
