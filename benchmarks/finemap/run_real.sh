#!/bin/bash
# Real GWAS loci (FinnGen ENDO, PCOS, fibroids, ANM; LD from 1000G EUR). Starts after the simulation benchmark.
cd /mnt/hdd2/loc-tran/susie_rs
until [ -f bench/logs/phase1b.done ]; do sleep 60; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
B=/mnt/hdd2/loc-tran/susie_rs/susie-rs/target/release/susiers
RB=/mnt/hdd2/loc-tran/susie_rs/renv/bin/Rscript
mkdir -p real/res real/logs
for t in PCOS ENDO FIB ANM; do
  $B batch --manifest real/$t/manifest_factor.tsv --out real/res/${t}_rs_auto --threads 1 --ld-mode auto > real/logs/${t}_rs_auto.log 2>&1
  $B batch --manifest real/$t/manifest.tsv --out real/res/${t}_rs_sym --threads 1 --ld-mode sym > real/logs/${t}_rs_sym.log 2>&1
  R_LIBS_USER= $RB bench/run_susier.R real/$t/manifest.tsv real/res/${t}_susier_ob > real/logs/${t}_susier_ob.log 2>&1
  /mnt/hdd2/loc-tran/miniconda3/bin/python bench/run_pysusie.py real/$t/manifest.tsv real/res/${t}_pysusie default > real/logs/${t}_pysusie.log 2>&1
done
# throughput on all 396 real loci, 48 pinned threads
cat real/*/manifest_factor.tsv > real/all_factor.tsv
$B batch --manifest real/all_factor.tsv --out real/res/all_rs_48 --threads 48 --pin --ld-mode auto > real/logs/all_rs_48.log 2>&1
echo DONE > real/real.done
