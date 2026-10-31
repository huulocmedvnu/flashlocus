# flashLocus

Summary-statistic fine-mapping and colocalisation in Rust. flashLocus is a line-by-line port of
susieR 0.14.2 (`susie_rss`) and coloc 5.2.3 (`coloc.bf_bf`, as used by `coloc.susie`), with exact
LD computation from PLINK 1 binary files.

## Install

    pip install flashlocus       # command-line tool + Python module (wheels: Linux x86_64/aarch64, macOS)

From source:

    cargo build --release        # binary: target/release/flashlocus
    cargo test --release         # regression tests against susieR, coloc and PLINK conventions
    pip install maturin && maturin develop --release   # Python module from a checkout

Requires Rust 1.75 or later to build. AVX-512 or AVX2 kernels are selected at run time (`flashlocus cpu`).

## Commands

| command | purpose |
|---|---|
| `rss` | SuSiE-RSS on one locus: z scores, LD matrix (p x p f64), sample size |
| `batch` | many loci from a manifest, parallel across loci, optional core pinning |
| `ld` | LD from `.bed/.bim/.fam` for a region or an ID list (`--extract`) |
| `prep` | GWAS windows to allele-aligned z scores and reference genotypes |
| `atlas` | coloc.susie between GWAS fits and eQTL Catalogue SuSiE lbf files |
| `sim`, `compare` | simulation and benchmarking helpers |

Examples:

    flashlocus rss --z locus.z --ld locus.ld.bin --n 150000 --out fit/locus
    flashlocus batch --manifest loci.tsv --out fits --threads 48 --pin --ld-mode auto
    flashlocus ld --bfile ref --chr 11 --start 29700000 --end 30700000 --extract ids.txt --out ld/fshb
    flashlocus atlas --loci loci.tsv \
        --qtl uterus=QTD000346.lbf_variable.txt.gz,QTD000346.credible_sets.tsv.gz --out atlas

Outputs per fit: `.pip`, `.cs.tsv` (effect, size, coverage, purity, members), `.lbf_variable.tsv`,
`.info.tsv`. `atlas` writes one `<tissue>.coloc.tsv` with PP.H0 to PP.H4 per credible-set pair.

## Python

```python
import numpy as np, flashlocus

fit = flashlocus.susie_rss(z, R, n=150000)            # z: (p,), R: (p, p); susieR 0.14.2 defaults
fit["pip"], fit["sets"][0]["members"]                  # members are 0-based
ids, chrom, pos, a1, a2, af, R = flashlocus.ld("ref", 11, 29_700_000, 30_700_000)
pairs = flashlocus.coloc_bf_bf(fit1["lbf_variable"][[0]], ids1, fit2["lbf_variable"][[0]], ids2)
pairs[0]["PP.H4"]
```

## Conventions

- SuSiE defaults follow susieR 0.14.2 (L = 10, coverage 0.95, min_abs_corr 0.5, prior variance
  optimised with Brent on log V, null check, ELBO tolerance 1e-3).
- coloc defaults follow coloc 5.2.3 (p1 = p2 = 1e-4, p12 = 5e-6, overlap.min 0.5), including its
  treatment of the last column in `logbf_to_pp` and the column order of the second dataset.
- LD: Pearson r of A1 dosages. Without missing calls the sums are exact integers (as in PLINK 1.9);
  with missing calls, missing dosages are mean-imputed. On chrX, male calls are haploid as in
  PLINK 1.9 (`--x-diploid` codes them 0/2).

## Validation summary (details in the accompanying manuscript)

- SuSiE: PIP r = 1.0000 against susieR 0.14.2 on 200 simulated and 396 real loci, maximum PIP
  difference 6e-8 (simulated) and 3e-5 (real).
- coloc: 27,159 credible-set pairs across 6 GTEx tissues, identical pair set and lead variants,
  maximum PP.H4 difference 5e-12.
- LD: 396 windows (1000 Genomes EUR), maximum difference 4.4e-16 against PLINK 1.9 and PLINK 2
  on autosomes and against PLINK 1.9 on chrX.

## License

MIT
