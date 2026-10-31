//! flashlocus: SuSiE-RSS fine-mapping (port of susieR 0.14.2), LD from PLINK, simulation and benchmarking.
//!   flashlocus rss     --z z.txt --ld ld.bin --n N [--L 10] [--f32] --out prefix
//!   flashlocus batch   --manifest loci.tsv --out dir [--threads T] [--pin] [--f32] [--L 10]
//!   flashlocus ld      --bfile ref --chr C --start S --end E [--maf 0.01] [--extract ids] [--x-diploid] [--float] --out prefix
//!                   chrX: male calls haploid as in PLINK 1.9 unless --x-diploid
//!   flashlocus sim     --bfile ref --out dir [--loci 200] [--seed 20261031]
//!   flashlocus compare --manifest loci.tsv --rust dir --ref dir
//!   flashlocus atlas   --loci loci.tsv --qtl name=lbf.txt.gz,credible_sets.tsv.gz [...] --out dir (coloc.susie)
//!   flashlocus prep    --bfile ref --windows list.tsv --format ssf|finngen|anm --match rsid|chrpos [--n N] --out dir
//! ld.bin: p*p little-endian f64, row-major. manifest: id<TAB>z_path<TAB>ld_path<TAB>n
use crate::ldmat::LdMat;
use crate::{atlas, bfile, compare, ldmat, prep, sim, simd, susie, topo};
use std::io::Write;
use std::time::Instant;

fn arg(a: &[String], k: &str) -> Option<String> { a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned()) }
fn flag(a: &[String], k: &str) -> bool { a.iter().any(|x| x == k) }
fn read_z(p: &str) -> Vec<f64> {
    std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{p}: {e}")).lines().filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse::<f64>().unwrap_or(f64::NAN)).collect()
}

fn write_fit(f: &susie::Fit, out: &str, secs: f64) {
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}.pip")).unwrap());
    for x in &f.pip { writeln!(w, "{:.10e}", x).unwrap() }
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}.cs.tsv")).unwrap());
    writeln!(w, "effect\tsize\tcoverage\tmin_abs_corr\tmean_abs_corr\tmedian_abs_corr\tmembers").unwrap();
    for c in &f.sets {
        writeln!(w, "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}", c.effect + 1, c.members.len(), c.coverage, c.purity[0], c.purity[1], c.purity[2],
                 c.members.iter().map(|x| (x + 1).to_string()).collect::<Vec<_>>().join(",")).unwrap();
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}.lbf_variable.tsv")).unwrap());
    for k in 0..f.l { writeln!(w, "{}", f.lbf_variable[k * f.p..(k + 1) * f.p].iter().map(|x| format!("{:.8e}", x)).collect::<Vec<_>>().join("\t")).unwrap() }
    let mut w = std::fs::File::create(format!("{out}.info.tsv")).unwrap();
    writeln!(w, "niter\t{}\nconverged\t{}\nelbo\t{:.8}\nseconds\t{:.6}\nV\t{}", f.niter, f.converged, f.elbo.last().unwrap_or(&f64::NAN), secs,
             f.v.iter().map(|x| format!("{:.6e}", x)).collect::<Vec<_>>().join(",")).unwrap();
}

fn fit_one<M: LdMat>(z: &[f64], ld: &M, n: f64, prm: &susie::Params) -> (Result<susie::Fit, String>, f64) {
    let t = Instant::now(); let f = susie::susie_rss(z, ld, n, prm); (f, t.elapsed().as_secs_f64())
}

/// LD storage modes: full (row-major f64), sym (packed upper f64, default), sym32 (packed upper f32)
fn run_mode(mode: &str, z: &[f64], lp: &str, n: f64, prm: &susie::Params) -> (Result<susie::Fit, String>, f64) {
    match mode {
        "full" => fit_one(z, &ldmat::read_full(lp, z.len()).unwrap(), n, prm),
        "sym32" => fit_one(z, &ldmat::read_sym_f32(lp, z.len()).unwrap(), n, prm),
        // auto: genotype input -> factor when p > 1.5 N, otherwise packed symmetric R built from it;
        // LD-matrix input -> packed symmetric f64
        "auto" => {
            if let Some((g, nn)) = lp.rsplit_once(':').filter(|(g, _)| g.ends_with(".geno.bin")) {
                let nn: usize = nn.parse().unwrap(); let f = ldmat::read_factor(g, z.len(), nn).unwrap();
                if z.len() as f64 > 1.5 * nn as f64 { fit_one(z, &f, n, prm) } else { fit_one(z, &ldmat::sym_from_factor(&f), n, prm) }
            } else { fit_one(z, &ldmat::read_sym_f64(lp, z.len()).unwrap(), n, prm) }
        }
        // factor modes: lp is <geno.bin>:<N> (standardised reference genotypes), R is never formed
        "factor" | "factor32" => {
            let (g, nn) = lp.rsplit_once(':').expect("factor mode needs geno.bin:N");
            let nn: usize = nn.parse().unwrap();
            if mode == "factor" { fit_one(z, &ldmat::read_factor(g, z.len(), nn).unwrap(), n, prm) } else { fit_one(z, &ldmat::read_factor_f32(g, z.len(), nn).unwrap(), n, prm) }
        }
        _ => fit_one(z, &ldmat::read_sym_f64(lp, z.len()).unwrap(), n, prm),
    }
}

/// Command-line entry point; `a[0]` is the program name.
pub fn run(a: Vec<String>) {
    let mut prm = susie::Params::default();
    if let Some(l) = arg(&a, "--L") { prm.l = l.parse().unwrap() }
    let mode = arg(&a, "--ld-mode").unwrap_or_else(|| if flag(&a, "--f32") { "sym32".into() } else { "auto".into() });
    match a.get(1).map(|s| s.as_str()) {
        Some("rss") => {
            let z = read_z(&arg(&a, "--z").expect("--z")); let lp = arg(&a, "--ld").expect("--ld");
            let n: f64 = arg(&a, "--n").expect("--n").parse().unwrap(); let out = arg(&a, "--out").expect("--out");
            let (f, s) = run_mode(&mode, &z, &lp, n, &prm);
            match f { Ok(f) => write_fit(&f, &out, s), Err(e) => { eprintln!("error: {e}"); std::process::exit(1) } }
        }
        Some("batch") => {
            use rayon::prelude::*;
            let man = std::fs::read_to_string(arg(&a, "--manifest").expect("--manifest")).unwrap();
            let out = arg(&a, "--out").expect("--out"); std::fs::create_dir_all(&out).unwrap();
            let threads: usize = arg(&a, "--threads").map(|x| x.parse().unwrap()).unwrap_or(1);
            let pool = topo::pool(threads, flag(&a, "--pin"));
            // relative paths in the manifest are resolved against the manifest's directory
            let mdir = std::path::Path::new(&arg(&a, "--manifest").unwrap()).parent().map(|p| p.to_path_buf()).unwrap_or_default();
            let abs = |s: &str| -> String { let (f, rest) = s.split_once(':').map(|(f, r)| (f, format!(":{r}"))).unwrap_or((s, String::new()));
                let p = std::path::Path::new(f); if p.is_absolute() { s.to_string() } else { format!("{}{rest}", mdir.join(p).display()) } };
            let rows: Vec<Vec<String>> = man.lines().filter(|l| !l.is_empty() && !l.starts_with('#')).map(|l| {
                let mut f: Vec<String> = l.split('\t').map(String::from).collect(); f[1] = abs(&f[1]); f[2] = abs(&f[2]); f }).collect();
            if let Some(r) = rows.first() { if !std::path::Path::new(&r[1]).exists() { eprintln!("error: {} not found", r[1]); std::process::exit(1) } }
            let t0 = Instant::now();
            let res: Vec<String> = pool.install(|| rows.par_iter().map(|r| {
                if !std::path::Path::new(&r[1]).exists() { return format!("{}\tNA\tNA\tNA\terror: missing {}", r[0], r[1]) }
                let z = read_z(&r[1]); let n: f64 = r[3].parse().unwrap();
                // LD is read inside the worker thread: first touch places it on the worker's NUMA node
                let (f, s) = run_mode(&mode, &z, &r[2], n, &prm);
                match f {
                    Ok(f) => { write_fit(&f, &format!("{out}/{}", r[0]), s); format!("{}\t{}\t{:.6}\t{}\t{}", r[0], z.len(), s, f.niter, f.sets.len()) }
                    Err(e) => format!("{}\t{}\tNA\tNA\terror: {}", r[0], z.len(), e),
                }
            }).collect());
            let mut w = std::fs::File::create(format!("{out}/timing.tsv")).unwrap();
            writeln!(w, "id\tp\tseconds\tniter\tn_cs").unwrap();
            for l in res { writeln!(w, "{l}").unwrap() }
            eprintln!("batch: {} loci in {:.2} s, {} threads, simd {}, ld {}", rows.len(), t0.elapsed().as_secs_f64(), threads, simd::level_name(), mode);
        }
        Some("ld") => {
            let b = bfile::Bfile::open(&arg(&a, "--bfile").expect("--bfile")).unwrap();
            let chr: u8 = arg(&a, "--chr").expect("--chr").parse().unwrap();
            let (s, e): (u32, u32) = (arg(&a, "--start").unwrap().parse().unwrap(), arg(&a, "--end").unwrap().parse().unwrap());
            let maf: f64 = arg(&a, "--maf").map(|x| x.parse().unwrap()).unwrap_or(0.01);
            let out = arg(&a, "--out").expect("--out");
            // --extract ids.txt restricts the window to listed variant IDs (bim order kept, MAF filter still applied)
            let ext: Option<std::collections::HashSet<String>> = arg(&a, "--extract").map(|f| std::fs::read_to_string(f).unwrap().lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect());
            let t = Instant::now();
            let (keep, af, r, mode) = match b.ld_region(chr, s, e, ext.as_ref(), maf, !flag(&a, "--x-diploid"), flag(&a, "--float")) {
                Some(x) => x, None => { eprintln!("no variants"); std::process::exit(1) }
            };
            let p = keep.len();
            let tc = t.elapsed().as_secs_f64();
            bfile::write_f64(&format!("{out}.ld.bin"), &r);
            let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{out}.snps.tsv")).unwrap());
            writeln!(w, "id\tchr\tpos\ta1\ta2\taf_a1").unwrap();
            for (k, &i) in keep.iter().enumerate() { let (id, a1, a2) = b.info(i); writeln!(w, "{id}\t{}\t{}\t{a1}\t{a2}\t{:.6}", b.chr[i], b.pos[i], af[k]).unwrap() }
            eprintln!("ld: {p} variants, {mode}, compute {tc:.3} s, total {:.2} s (simd {})", t.elapsed().as_secs_f64(), simd::level_name());
        }
        Some("sim") => {
            let loci: usize = arg(&a, "--loci").map(|x| x.parse().unwrap()).unwrap_or(200);
            let seed: u64 = arg(&a, "--seed").map(|x| x.parse().unwrap()).unwrap_or(20261031);
            sim::run(&arg(&a, "--bfile").expect("--bfile"), &arg(&a, "--out").expect("--out"), loci, seed, &[500, 1000, 2000, 3000, 5000]);
        }
        Some("compare") => compare::run(&arg(&a, "--manifest").expect("--manifest"), &arg(&a, "--rust").expect("--rust"), &arg(&a, "--ref").expect("--ref")),
        Some("prep") => prep::run(&arg(&a, "--bfile").expect("--bfile"), &arg(&a, "--windows").expect("--windows"), &arg(&a, "--format").unwrap_or("ssf".into()),
                                 arg(&a, "--match").map(|x| x == "rsid").unwrap_or(true), arg(&a, "--n").map(|x| x.parse().unwrap()),
                                 arg(&a, "--maf").map(|x| x.parse().unwrap()).unwrap_or(0.01), &arg(&a, "--out").expect("--out")),
        Some("atlas") => atlas::run(&a[2..]),
        Some("cpu") => { println!("simd\t{}\nphysical_cores\t{:?}", simd::level_name(), topo::physical_cores()) }
        _ => { eprintln!("usage: flashlocus rss|batch|ld|sim|compare|prep|atlas|cpu ..."); std::process::exit(2) }
    }
}
