#!/bin/bash
cd ${BENCH_ROOT:-.}/ldbench/chrx
R=${REF_1KG}; B=../../susie-rs/target/release/susiers; P1=${PLINK19:-plink}; PY=python
for id in ANM0236 ANM0237 ANM0238 ANM0239 ANM0240 ANM0241 ANM0242 ANM0001 PCOS0001 ENDO0001; do
  t=${id//[0-9]/}; s=../../real/$t/$id.snps.tsv; [ -f $s ] || continue
  tail -n +2 $s | cut -f1 > w.ids; p=$(wc -l < w.ids); c=$(sed -n 2p $s | cut -f2)
  lo=$(tail -n +2 $s | cut -f3 | sort -n | head -1); hi=$(tail -n +2 $s | cut -f3 | sort -n | tail -1)
  $B ld --bfile $R --chr $c --start $lo --end $hi --maf 0 --extract w.ids --out e 2> e.err
  $B ld --bfile $R --chr $c --start $lo --end $hi --maf 0 --extract w.ids --float --out f 2> f.err
  $P1 --bfile $R --extract w.ids --keep-allele-order --r square bin --threads 48 --out q > /dev/null 2>&1
  $PY -I -c "
import numpy as np
a=np.fromfile('e.ld.bin'); f=np.fromfile('f.ld.bin'); q=np.fromfile('q.ld.bin')
print('$id chr$c p=$p', '$(grep -o 'exact-int\|float' e.err)', 'exact vs plink1.9 %.2e'%np.nanmax(np.abs(a-q)), 'float vs plink1.9 %.2e'%np.nanmax(np.abs(f-q)))"
done
