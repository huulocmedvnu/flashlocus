#!/bin/bash
# all tools single-threaded
cd ${BENCH_ROOT:-.}/bench
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
M=sim/manifest.tsv
${BENCH_ROOT:-.}/susie-rs/target/release/susiers batch --manifest $M --out res/rust --threads 1 > logs/rust.log 2>&1
Rscript run_susier.R $M res/susier > logs/susier.log 2>&1
python run_pysusie.py $M res/pysusie_default default > logs/py_default.log 2>&1
python run_pysusie.py $M res/pysusie_exact exact > logs/py_exact.log 2>&1
echo PHASE1_DONE > logs/phase1.done
