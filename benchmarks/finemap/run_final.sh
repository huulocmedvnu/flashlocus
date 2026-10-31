#!/bin/bash
# Final benchmark under the fairness rules in PLAN.md. Single-thread runs pinned to core 2; throughput on 48 cores.
cd /mnt/hdd2/loc-tran/susie_rs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
B=/mnt/hdd2/loc-tran/susie_rs/susie-rs/target/release/susiers
RB=/mnt/hdd2/loc-tran/susie_rs/renv/bin/Rscript
PY=/mnt/hdd2/loc-tran/miniconda3/bin/python
L016=/mnt/hdd2/loc-tran/susie_rs/renv_lib016
mkdir -p final/logs
run_set () {   # $1 name, $2 ld manifest, $3 factor manifest
  n=$1; M=$2; MF=$3; o=final/$n
  # preflight: every file of the first locus must exist, else stop (a run on missing files once produced empty timings)
  for f in "$(head -1 $M | cut -f2)" "$(head -1 $M | cut -f3)" "$(head -1 $MF | cut -f3 | cut -d: -f1)"; do
    [ -f "$f" ] || { echo "PREFLIGHT FAIL: $f" > final/logs/${n}_preflight.log; exit 1; }
  done
  taskset -c 2 $B batch --manifest $M --out $o/rs_sym --threads 1 --ld-mode sym > final/logs/${n}_rs_sym.log 2>&1
  taskset -c 2 $B batch --manifest $M --out $o/rs_sym32 --threads 1 --ld-mode sym32 > final/logs/${n}_rs_sym32.log 2>&1
  taskset -c 2 $B batch --manifest $MF --out $o/rs_auto --threads 1 --ld-mode auto > final/logs/${n}_rs_auto.log 2>&1
  taskset -c 2 $RB bench/run_susier.R $M $o/susier014 > final/logs/${n}_susier014.log 2>&1
  SUSIE_LIB=$L016 taskset -c 2 $RB bench/run_susier.R $M $o/susier016 > final/logs/${n}_susier016.log 2>&1
  taskset -c 2 $PY bench/run_pysusie.py $M $o/pysusie default > final/logs/${n}_pysusie.log 2>&1
  taskset -c 2 $PY bench/run_pysusie.py $M $o/pysusie_exact exact > final/logs/${n}_pysusie_exact.log 2>&1
  # throughput, 48 physical cores
  { /usr/bin/time -f "susiers_sym wall %e" $B batch --manifest $M --out $o/tp_rs_sym --threads 48 --pin --ld-mode sym ; } > final/logs/${n}_tp.log 2>&1
  { /usr/bin/time -f "susiers_auto wall %e" $B batch --manifest $MF --out $o/tp_rs_auto --threads 48 --pin --ld-mode auto ; } >> final/logs/${n}_tp.log 2>&1
  $RB bench/run_susier_par.R $M 48 >> final/logs/${n}_tp.log 2>&1
  SUSIE_LIB=$L016 $RB bench/run_susier_par.R $M 48 >> final/logs/${n}_tp.log 2>&1
  $PY bench/run_pysusie_par.py $M 48 >> final/logs/${n}_tp.log 2>&1
}
run_set sim bench/sim2/manifest.tsv bench/sim2/manifest_factor.tsv
cat real/{PCOS,ENDO,FIB,ANM}/manifest.tsv > real/all.tsv
cat real/{PCOS,ENDO,FIB,ANM}/manifest_factor.tsv > real/all_factor.tsv
run_set real real/all.tsv real/all_factor.tsv
echo DONE > final/final.done
