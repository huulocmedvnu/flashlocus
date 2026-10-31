//! coloc.susie atlas: GWAS flashlocus fits x eQTL Catalogue SuSiE lbf files.
//!   flashlocus atlas --loci loci.tsv --qtl name=lbf.txt.gz,credible_sets.tsv.gz [--qtl ...] --out dir
//!                 [--threads T] [--p1 1e-4] [--p2 1e-4] [--p12 5e-6] [--overlap-min 0.5]
//! loci.tsv: id<TAB>fit_prefix<TAB>keys_path. keys file: one line per GWAS SNP (fit order),
//! chr_pos38_alleleA_alleleB with alleles sorted, or NA (never matches).
//! A gene is tested against a locus when its cis region overlaps [min, max] GRCh38 position of the locus.
use crate::coloc;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};
use std::time::Instant;

struct Locus { id: String, chr: String, lo: u64, hi: u64, keys: Vec<String>, effects: Vec<usize>, rows: Vec<Vec<f64>> }

struct Gene { id: String, chr: String, lo: u64, hi: u64, ls: Vec<usize>, keys: Vec<String>, rows: Vec<Vec<f64>> }

fn read_locus(id: &str, prefix: &str, keys_path: &str) -> Option<Locus> {
    let keys: Vec<String> = std::fs::read_to_string(keys_path).unwrap_or_else(|e| panic!("{keys_path}: {e}"))
        .lines().enumerate().map(|(i, l)| if l == "NA" { format!("NA#{i}") } else { l.to_string() }).collect();
    let cs = std::fs::read_to_string(format!("{prefix}.cs.tsv")).unwrap_or_else(|e| panic!("{prefix}.cs.tsv: {e}"));
    let effects: Vec<usize> = cs.lines().skip(1).filter(|l| !l.is_empty()).map(|l| l.split('\t').next().unwrap().parse::<usize>().unwrap() - 1).collect();
    if effects.is_empty() { return None }
    let lbf = std::fs::read_to_string(format!("{prefix}.lbf_variable.tsv")).unwrap();
    let all: Vec<&str> = lbf.lines().collect();
    let rows: Vec<Vec<f64>> = effects.iter().map(|&e| all[e].split('\t').map(|x| x.parse::<f64>().unwrap()).collect()).collect();
    assert!(rows.iter().all(|r| r.len() == keys.len()), "{id}: keys and lbf width differ");
    let (mut chr, mut lo, mut hi) = (String::new(), u64::MAX, 0u64);
    for k in keys.iter().filter(|k| !k.starts_with("NA#")) {
        let mut f = k.split('_');
        let c = f.next().unwrap(); let p: u64 = f.next().unwrap().parse().unwrap();
        if chr.is_empty() { chr = c.to_string() }
        if c == chr { lo = lo.min(p); hi = hi.max(p) }
    }
    if chr.is_empty() { return None }
    Some(Locus { id: id.to_string(), chr, lo, hi, keys, effects, rows })
}

fn parse_region(r: &str) -> (String, u64, u64) {
    let (c, ab) = r.split_once(':').unwrap();
    let (a, b) = ab.rsplit_once('-').unwrap();  // start can be negative near chromosome ends ("chr1:-80000-1920000")
    (c.trim_start_matches("chr").to_string(), a.parse::<i64>().unwrap().max(0) as u64, b.parse().unwrap())
}

/// credible_sets.tsv.gz -> gene -> (region, sorted CS L indices, 0-based)
fn read_cs(path: &str) -> HashMap<String, (String, Vec<usize>)> {
    let r = BufReader::new(flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"))));
    let mut lines = r.lines();
    let h: Vec<String> = lines.next().unwrap().unwrap().split('\t').map(|s| s.to_string()).collect();
    let ix = |c: &str| h.iter().position(|x| x == c).unwrap_or_else(|| panic!("{path}: no column {c}"));
    let (it, ic, ir) = (ix("molecular_trait_id"), ix("cs_id"), ix("region"));
    let mut m: HashMap<String, (String, Vec<usize>)> = HashMap::new();
    for l in lines {
        let l = l.unwrap(); let f: Vec<&str> = l.split('\t').collect();
        let k: usize = f[ic].rsplit_once("_L").unwrap().1.parse().unwrap();
        let e = m.entry(f[it].to_string()).or_insert_with(|| (f[ir].to_string(), vec![]));
        if !e.1.contains(&(k - 1)) { e.1.push(k - 1) }
    }
    for v in m.values_mut() { v.1.sort() }
    m
}

fn variant_key(v: &str) -> String {
    let mut f = v.split('_');
    let c = f.next().unwrap().trim_start_matches("chr"); let p = f.next().unwrap();
    let (a, b) = (f.next().unwrap(), f.next().unwrap());
    if a <= b { format!("{c}_{p}_{a}_{b}") } else { format!("{c}_{p}_{b}_{a}") }
}

type LbfRow = (usize, String, Vec<f64>);

fn parse_chunk(buf: &[u8], genes: &HashMap<String, usize>, ls: &[Vec<usize>]) -> Vec<LbfRow> {
    let mut out = Vec::new();
    for line in buf.split(|&b| b == b'\n') {
        if line.is_empty() { continue }
        let tab = match line.iter().position(|&b| b == b'\t') { Some(t) => t, None => continue };
        let g = match std::str::from_utf8(&line[..tab]).ok().and_then(|t| genes.get(t)) { Some(&g) => g, None => continue };
        let s = std::str::from_utf8(line).unwrap();
        let f: Vec<&str> = s.split('\t').collect();
        if f[0] == "molecular_trait_id" { continue }
        let v: Vec<f64> = ls[g].iter().map(|&k| f[5 + k].parse::<f64>().unwrap()).collect();
        out.push((g, variant_key(f[2]), v));
    }
    out
}

fn check_header(h: &[u8], path: &str) {
    let h = std::str::from_utf8(h).unwrap();
    let hc: Vec<&str> = h.trim_end().split('\t').collect();
    assert!(hc.len() > 5 && hc[0] == "molecular_trait_id" && hc[2] == "variant" && hc[5] == "lbf_variable1", "{path}: unexpected header");
}

/// Default: whole gzip member decompressed in memory with libdeflate (about twice zlib speed), then lines parsed
/// on the rayon pool in file order. FLASHLOCUS_STREAM=1: streaming flate2 (zlib-rs) with bounded memory.
fn read_lbf(path: &str, genes: &HashMap<String, usize>, ls: &[Vec<usize>]) -> Vec<LbfRow> {
    if std::env::var("FLASHLOCUS_STREAM").is_ok() { return read_lbf_stream(path, genes, ls) }
    let prof = std::env::var("FLASHLOCUS_PROFILE").is_ok(); let t = Instant::now();
    let gz = std::fs::read(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    if prof { eprintln!("  read gz {:.2} s", t.elapsed().as_secs_f64()) }
    let isize = u32::from_le_bytes(gz[gz.len() - 4..].try_into().unwrap()) as usize;  // output size mod 2^32
    // gzip ISIZE is the output size mod 2^32. Start from the smallest candidate >= 2.5x the compressed size
    // (text compresses at least that well); a larger buffer is harmless because the true length is returned.
    let mut d = libdeflater::Decompressor::new();
    let mut k = 0usize;
    while isize + (k << 32) < gz.len() * 5 / 2 { k += 1 }
    let data = loop {
        let cap = isize + (k << 32);
        let mut buf: Vec<u8> = Vec::with_capacity(cap);
        // SAFETY: libdeflate writes the output region and returns its length; only that prefix is kept.
        unsafe { buf.set_len(cap) }
        match d.gzip_decompress(&gz, &mut buf) {
            Ok(n) => { buf.truncate(n); break buf }
            Err(libdeflater::DecompressionError::InsufficientSpace) => { k += 1; assert!(k < 64, "{path}: cannot size output") }
            Err(e) => panic!("{path}: {e:?}"),
        }
    };
    drop(gz);
    if prof { eprintln!("  inflate {:.2} s ({} bytes)", t.elapsed().as_secs_f64(), data.len()) }
    let hend = data.iter().position(|&b| b == b'\n').unwrap();
    check_header(&data[..hend], path);
    let mut cuts = vec![hend + 1];
    let step = 1 << 23;
    while *cuts.last().unwrap() < data.len() {
        let a = *cuts.last().unwrap();
        let mut b = (a + step).min(data.len());
        while b < data.len() && data[b - 1] != b'\n' { b += 1 }
        cuts.push(b);
    }
    let parts: Vec<Vec<LbfRow>> = cuts.par_windows(2).map(|w| parse_chunk(&data[w[0]..w[1]], genes, ls)).collect();
    if prof { eprintln!("  parse {:.2} s", t.elapsed().as_secs_f64()) }
    let r: Vec<LbfRow> = parts.into_iter().flatten().collect();
    if prof { eprintln!("  flatten {:.2} s", t.elapsed().as_secs_f64()) }
    r
}

fn read_lbf_stream(path: &str, genes: &HashMap<String, usize>, ls: &[Vec<usize>]) -> Vec<LbfRow> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut r = BufReader::with_capacity(1 << 22, flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 22, f)));
    let mut h = Vec::new(); r.read_until(b'\n', &mut h).unwrap();
    check_header(&h, path);
    let (tx, rx) = std::sync::mpsc::sync_channel::<(usize, Vec<u8>)>(96);
    let reader = std::thread::spawn(move || {
        let mut seq = 0;
        loop {
            let mut buf = vec![0u8; 1 << 23];
            let mut n = 0;
            while n < buf.len() { let k = r.read(&mut buf[n..]).unwrap(); if k == 0 { break } n += k }
            buf.truncate(n);
            if n == 0 { break }
            r.read_until(b'\n', &mut buf).unwrap();
            tx.send((seq, buf)).unwrap(); seq += 1;
        }
    });
    let mut parts: Vec<(usize, Vec<LbfRow>)> = rx.into_iter().par_bridge().map(|(seq, buf)| (seq, parse_chunk(&buf, genes, ls))).collect();
    reader.join().unwrap();
    parts.sort_by_key(|p| p.0);
    parts.into_iter().flat_map(|p| p.1).collect()
}

fn load_qtl(lbf: &str, cs: &str, loci: &[Locus]) -> Vec<Gene> {
    let t = Instant::now();
    let csm = read_cs(cs);
    if std::env::var("FLASHLOCUS_PROFILE").is_ok() { eprintln!("  credible sets {:.2} s", t.elapsed().as_secs_f64()) }
    let mut genes: Vec<Gene> = csm.into_iter().filter_map(|(id, (reg, ls))| {
        let (chr, lo, hi) = parse_region(&reg);
        if loci.iter().any(|l| l.chr == chr && lo <= l.hi && hi >= l.lo) { Some(Gene { id, chr, lo, hi, ls, keys: vec![], rows: vec![] }) } else { None }
    }).collect();
    genes.sort_by(|a, b| a.id.cmp(&b.id));
    let idx: HashMap<String, usize> = genes.iter().enumerate().map(|(i, g)| (g.id.clone(), i)).collect();
    let ls: Vec<Vec<usize>> = genes.iter().map(|g| g.ls.clone()).collect();
    let mut cols: Vec<Vec<Vec<f64>>> = vec![vec![]; genes.len()];
    for (g, key, v) in read_lbf(lbf, &idx, &ls) { genes[g].keys.push(key); cols[g].push(v) }
    for (g, c) in genes.iter_mut().zip(cols) {
        // duplicate variant keys within a gene: first kept, later ones renamed so they never match
        let mut seen = std::collections::HashSet::new();
        for (k, key) in g.keys.iter_mut().enumerate() { if !seen.insert(key.clone()) { *key = format!("DUP#{k}") } }
        g.rows = (0..g.ls.len()).map(|r| c.iter().map(|v| v[r]).collect()).collect();
    }
    genes.retain(|g| !g.keys.is_empty());
    genes
}

fn arg(a: &[String], k: &str) -> Option<String> { a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned()) }

pub fn run(a: &[String]) {
    let out = arg(a, "--out").expect("--out");
    std::fs::create_dir_all(&out).unwrap();
    let pf = |k: &str, d: f64| arg(a, k).map(|v| v.parse::<f64>().unwrap()).unwrap_or(d);
    let (p1, p2, p12, omin) = (pf("--p1", 1e-4), pf("--p2", 1e-4), pf("--p12", 5e-6), pf("--overlap-min", 0.5));
    let t0 = Instant::now();
    let list = std::fs::read_to_string(arg(a, "--loci").expect("--loci")).unwrap();
    let rows: Vec<Vec<String>> = list.lines().filter(|l| !l.is_empty()).map(|l| l.split('\t').map(|s| s.to_string()).collect()).collect();
    let loci: Vec<Locus> = rows.par_iter().filter_map(|r| read_locus(&r[0], &r[1], &r[2])).collect();
    eprintln!("atlas: {} loci with a credible set ({:.1} s)", loci.len(), t0.elapsed().as_secs_f64());
    let qtls: Vec<(String, String, String)> = a.iter().enumerate().filter(|(_, x)| *x == "--qtl").map(|(i, _)| {
        let (n, p) = a[i + 1].split_once('=').unwrap(); let (l, c) = p.split_once(',').unwrap();
        (n.to_string(), l.to_string(), c.to_string())
    }).collect();
    let mut log = std::fs::File::create(format!("{out}/timing.tsv")).unwrap();
    writeln!(log, "tissue\tgenes\tpairs\tunaligned_bf2\tload_s\tcoloc_s").unwrap();
    // Tissues are loaded concurrently (each gzip member inflates on one thread; parsing shares the pool).
    // FLASHLOCUS_QTL_SEQ=1 loads them one at a time.
    let seq = std::env::var("FLASHLOCUS_QTL_SEQ").is_ok();
    let tl0 = Instant::now();
    let loaded: Vec<(Vec<Gene>, f64)> = if seq {
        qtls.iter().map(|(_, l, c)| { let t = Instant::now(); let g = load_qtl(l, c, &loci); (g, t.elapsed().as_secs_f64()) }).collect()
    } else {
        qtls.par_iter().map(|(_, l, c)| { let t = Instant::now(); let g = load_qtl(l, c, &loci); (g, t.elapsed().as_secs_f64()) }).collect()
    };
    eprintln!("atlas: {} tissues loaded in {:.1} s ({})", qtls.len(), tl0.elapsed().as_secs_f64(), if seq { "sequential" } else { "concurrent" });
    for ((name, _, _), (genes, tl)) in qtls.iter().zip(loaded) {
        let tl = tl;
        let t = Instant::now();
        let res: Vec<(String, usize)> = loci.par_iter().map(|l| {
            let mut s = String::new(); let mut un = 0;
            let b1: Vec<&[f64]> = l.rows.iter().map(|r| r.as_slice()).collect();
            for g in genes.iter().filter(|g| g.chr == l.chr && g.lo <= l.hi && g.hi >= l.lo) {
                let b2: Vec<&[f64]> = g.rows.iter().map(|r| r.as_slice()).collect();
                let o = coloc::bf_bf(&b1, &l.keys, &b2, &g.keys, p1, p2, p12, omin);
                if o.bf2_unaligned { un += 1 }
                for p in &o.pairs {
                    s += &format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\n", l.id, name, g.id,
                                  l.effects[p.i] + 1, g.ls[p.j] + 1, p.nsnps, l.keys[p.hit1], g.keys[p.hit2], p.pp[0], p.pp[1], p.pp[2], p.pp[3], p.pp[4]);
                }
            }
            (s, un)
        }).collect();
        let tc = t.elapsed().as_secs_f64();
        let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}/{name}.coloc.tsv")).unwrap());
        writeln!(w, "locus\ttissue\tgene\tgwas_L\tqtl_L\tnsnps\thit1\thit2\tPP.H0\tPP.H1\tPP.H2\tPP.H3\tPP.H4").unwrap();
        let mut un = 0;
        for (s, u) in &res { w.write_all(s.as_bytes()).unwrap(); un += u }
        let np: usize = res.iter().map(|r| r.0.lines().count()).sum();
        writeln!(log, "{name}\t{}\t{np}\t{un}\t{tl:.3}\t{tc:.3}", genes.len()).unwrap();
        eprintln!("atlas {name}: {} genes, {np} pairs, load {tl:.1} s, coloc {tc:.2} s", genes.len());
    }
    eprintln!("atlas total {:.1} s", t0.elapsed().as_secs_f64());
}
