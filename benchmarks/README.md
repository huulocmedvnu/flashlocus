# Benchmarks

Scripts used for the validation and timing reported in the manuscript, run on a 2 x Xeon Platinum 8168
server (48 cores, AVX-512). They were run before the rename, so they call the binary `susiers`
(`../susie-rs/target/release/susiers`); the current binary is `flashlocus` with the same subcommands.

- `finemap/`: phase 1, SuSiE-RSS parity and speed against susieR 0.14.2 / 0.16.6 and pysusie
  (simulated 1000 Genomes EUR loci and 396 real loci from four reproductive GWAS). `run_final.sh` is the
  final fair benchmark (OpenBLAS + Rfast for R, every tool pinned to one core, then 48-core throughput).
- `coloc_atlas/`: phase 2, coloc.susie parity and speed against coloc 5.2.3 across six GTEx tissues from
  eQTL Catalogue (`run_atlas.sh`, R reference `atlas_r.R`, per-pair comparison `parity.py`).
- `ld/`: LD against PLINK 1.9 and PLINK 2 on the 396 windows (`ld_bench.sh`) and the chrX check.
