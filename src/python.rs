//! Python bindings (feature "python", built by maturin as flashlocus._core).
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::{bfile, coloc, ldmat, susie};

/// SuSiE-RSS fit (susieR 0.14.2 susie_rss defaults). z: (p,), R: (p, p) correlation, n: GWAS sample size.
/// Returns a dict with pip, alpha, lbf_variable, V, elbo, niter, converged and sets (0-based members).
#[pyfunction]
#[pyo3(signature = (z, R, n, L=10, coverage=0.95, min_abs_corr=0.5, max_iter=100, tol=1e-3,
                    scaled_prior_variance=0.2, estimate_prior_variance=true))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn susie_rss<'py>(py: Python<'py>, z: PyReadonlyArray1<'py, f64>, R: PyReadonlyArray2<'py, f64>, n: f64, L: usize,
                  coverage: f64, min_abs_corr: f64, max_iter: usize, tol: f64, scaled_prior_variance: f64,
                  estimate_prior_variance: bool) -> PyResult<Bound<'py, PyDict>> {
    let z = z.as_array().to_vec();
    let p = z.len();
    let ra = R.as_array();
    if ra.shape() != [p, p] { return Err(PyValueError::new_err(format!("R must be {p} x {p}, got {:?}", ra.shape()))) }
    let ld = ldmat::LdF64 { p, r: ra.iter().cloned().collect() };
    let prm = susie::Params { l: L, coverage, min_abs_corr, max_iter, tol, scaled_prior_variance, estimate_prior_variance, ..Default::default() };
    let f = py.detach(|| susie::susie_rss(&z, &ld, n, &prm)).map_err(PyRuntimeError::new_err)?;
    let d = PyDict::new(py);
    d.set_item("pip", f.pip.clone().into_pyarray(py))?;
    d.set_item("alpha", Array2::from_shape_vec((f.l, f.p), f.alpha.clone()).unwrap().into_pyarray(py))?;
    d.set_item("lbf_variable", Array2::from_shape_vec((f.l, f.p), f.lbf_variable.clone()).unwrap().into_pyarray(py))?;
    d.set_item("V", f.v.clone().into_pyarray(py))?;
    d.set_item("elbo", f.elbo.clone().into_pyarray(py))?;
    d.set_item("niter", f.niter)?;
    d.set_item("converged", f.converged)?;
    let sets = PyList::empty(py);
    for c in &f.sets {
        let s = PyDict::new(py);
        s.set_item("effect", c.effect)?;
        s.set_item("members", c.members.clone())?;
        s.set_item("coverage", c.coverage)?;
        s.set_item("min_abs_corr", c.purity[0])?;
        s.set_item("mean_abs_corr", c.purity[1])?;
        s.set_item("median_abs_corr", c.purity[2])?;
        sets.append(s)?;
    }
    d.set_item("sets", sets)?;
    Ok(d)
}

/// coloc.bf_bf (coloc 5.2.3): bf1 (k1, p1) and bf2 (k2, p2) log Bayes factors with SNP names snps1/snps2.
/// Returns one dict per tested row pair: idx1, idx2 (0-based rows), nsnps, hit1, hit2, PP.H0 ... PP.H4.
#[pyfunction]
#[pyo3(signature = (bf1, snps1, bf2, snps2, p1=1e-4, p2=1e-4, p12=5e-6, overlap_min=0.5))]
#[allow(clippy::too_many_arguments)]
fn coloc_bf_bf<'py>(py: Python<'py>, bf1: PyReadonlyArray2<'py, f64>, snps1: Vec<String>, bf2: PyReadonlyArray2<'py, f64>,
                    snps2: Vec<String>, p1: f64, p2: f64, p12: f64, overlap_min: f64) -> PyResult<Bound<'py, PyList>> {
    let rows = |a: &PyReadonlyArray2<'py, f64>, k: usize, name: &str| -> PyResult<Vec<Vec<f64>>> {
        let a = a.as_array();
        if a.ncols() != k { return Err(PyValueError::new_err(format!("{name} has {} columns but {k} SNP names", a.ncols()))) }
        Ok(a.rows().into_iter().map(|r| r.to_vec()).collect())
    };
    let (m1, m2) = (rows(&bf1, snps1.len(), "bf1")?, rows(&bf2, snps2.len(), "bf2")?);
    let r1: Vec<&[f64]> = m1.iter().map(|r| r.as_slice()).collect();
    let r2: Vec<&[f64]> = m2.iter().map(|r| r.as_slice()).collect();
    let o = coloc::bf_bf(&r1, &snps1, &r2, &snps2, p1, p2, p12, overlap_min);
    let out = PyList::empty(py);
    for p in &o.pairs {
        let d = PyDict::new(py);
        d.set_item("idx1", p.i)?; d.set_item("idx2", p.j)?; d.set_item("nsnps", p.nsnps)?;
        d.set_item("hit1", &snps1[p.hit1])?; d.set_item("hit2", &snps2[p.hit2])?;
        for h in 0..5 { d.set_item(format!("PP.H{h}"), p.pp[h])?; }
        out.append(d)?;
    }
    Ok(out)
}

/// LD (Pearson r of A1 dosages) from a PLINK 1 fileset for chrom:start-end, optionally restricted to IDs.
/// Exact integer sums when no call is missing; chrX males haploid as in PLINK 1.9 unless x_haploid=False.
/// Returns (ids, chrom, pos, a1, a2, af_a1, R).
#[pyfunction]
#[pyo3(signature = (bfile, chrom, start, end, extract=None, maf=0.0, x_haploid=true))]
#[allow(clippy::type_complexity)]
fn ld<'py>(py: Python<'py>, bfile: &str, chrom: u8, start: u32, end: u32, extract: Option<Vec<String>>, maf: f64, x_haploid: bool)
    -> PyResult<(Vec<String>, u8, Vec<u32>, Vec<String>, Vec<String>, Vec<f64>, Bound<'py, numpy::PyArray2<f64>>)> {
    let b = bfile::Bfile::open(bfile).map_err(PyRuntimeError::new_err)?;
    let ext: Option<std::collections::HashSet<String>> = extract.map(|v| v.into_iter().collect());
    let res = py.detach(|| b.ld_region(chrom, start, end, ext.as_ref(), maf, x_haploid, false));
    let (keep, af, r, _) = res.ok_or_else(|| PyValueError::new_err("no variants in region"))?;
    let p = keep.len();
    let (mut ids, mut a1s, mut a2s) = (vec![], vec![], vec![]);
    for &i in &keep { let (id, a1, a2) = b.info(i); ids.push(id); a1s.push(a1); a2s.push(a2) }
    let pos = keep.iter().map(|&i| b.pos[i]).collect();
    Ok((ids, chrom, pos, a1s, a2s, af, Array2::from_shape_vec((p, p), r).unwrap().into_pyarray(py)))
}

/// Run the command-line interface with argv (argv[0] is the program name).
#[pyfunction]
fn cli(argv: Vec<String>) { crate::cli::run(argv) }

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(susie_rss, m)?)?;
    m.add_function(wrap_pyfunction!(coloc_bf_bf, m)?)?;
    m.add_function(wrap_pyfunction!(ld, m)?)?;
    m.add_function(wrap_pyfunction!(cli, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
