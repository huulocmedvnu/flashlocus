#!/bin/bash
# all tools single-threaded
cd /mnt/hdd2/loc-tran/susie_rs/bench
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
M=sim/manifest.tsv
/mnt/hdd2/loc-tran/susie_rs/susie-rs/target/release/susiers batch --manifest $M --out res/rust --threads 1 > logs/rust.log 2>&1
Rscript run_susier.R $M res/susier > logs/susier.log 2>&1
/mnt/hdd2/loc-tran/miniconda3/bin/python run_pysusie.py $M res/pysusie_default default > logs/py_default.log 2>&1
/mnt/hdd2/loc-tran/miniconda3/bin/python run_pysusie.py $M res/pysusie_exact exact > logs/py_exact.log 2>&1
echo PHASE1_DONE > logs/phase1.done
