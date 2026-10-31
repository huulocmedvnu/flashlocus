#!/bin/bash
# Phase 1 on Rust-simulated loci (sim2): every tool single-threaded, then susiers thread scaling.
cd /mnt/hdd2/loc-tran/susie_rs/bench
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
B=/mnt/hdd2/loc-tran/susie_rs/susie-rs/target/release/susiers
M=sim2/manifest.tsv; MF=sim2/manifest_factor.tsv
for m in full sym sym32; do $B batch --manifest $M --out res2/rs_$m --threads 1 --ld-mode $m > logs/b_$m.log 2>&1; done
for m in factor auto; do $B batch --manifest $MF --out res2/rs_$m --threads 1 --ld-mode $m > logs/b_$m.log 2>&1; done
Rscript run_susier.R $M res2/susier > logs/b_susier.log 2>&1
/mnt/hdd2/loc-tran/miniconda3/bin/python run_pysusie.py $M res2/pysusie_default default > logs/b_pyd.log 2>&1
/mnt/hdd2/loc-tran/miniconda3/bin/python run_pysusie.py $M res2/pysusie_exact exact > logs/b_pye.log 2>&1
# throughput scaling (auto mode, genotype factors), pinned to physical cores vs not
for t in 1 6 12 24 48 96; do
  for pin in pin nopin; do
    f=""; [ $pin = pin ] && f="--pin"
    $B batch --manifest $MF --out res2/scale_${t}_$pin --threads $t $f --ld-mode auto > logs/scale_${t}_$pin.log 2>&1
  done
done
echo DONE > logs/phase1b.done
