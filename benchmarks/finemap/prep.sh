#!/bin/bash
cd ${BENCH_ROOT:-.}/real
export RAYON_NUM_THREADS=6
B=${BENCH_ROOT:-.}/susie-rs/target/release/susiers
R=${REF_1KG}
$B prep --bfile $R --windows win_PCOS.tsv --format ssf --match rsid --n 444734 --out PCOS > prep_PCOS.log 2>&1
$B prep --bfile $R --windows win_FIB.tsv --format ssf --match rsid --n 434152 --out FIB > prep_FIB.log 2>&1
$B prep --bfile $R --windows win_ENDO.tsv --format finngen --match rsid --n 150350 --out ENDO > prep_ENDO.log 2>&1
$B prep --bfile $R --windows win_ANM.tsv --format anm --match chrpos --out ANM > prep_ANM.log 2>&1
echo PREP_DONE > prep.done
