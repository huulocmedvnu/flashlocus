#!/bin/bash
# LD benchmark: 396 real windows, same SNP list per window, one invocation per window per tool, 48 threads each,
# reference files warm in page cache, wall time includes reading .bim/.bed and writing the p x p binary matrix.
cd ${BENCH_ROOT:-.}; W=ldbench/full; mkdir -p $W; cd $W
R=${REF_1KG}
B=../../susie-rs/target/release/susiers; P1=${PLINK19:-plink}; P2=../../tools/plink2/plink2; PY=python
cat $R.bed $R.bim $R.fam > /dev/null
echo -e "id\tp\tsusiers_s\tplink19_s\tplink2_s\tmaxdiff_p19\tmaxdiff_p2_abs" > ../ld_times.tsv
for s in ../../real/*/*.snps.tsv; do
  id=$(basename $s .snps.tsv)
  tail -n +2 $s | cut -f1 > w.ids; p=$(wc -l < w.ids); c=$(sed -n 2p $s | cut -f2)
  lo=$(tail -n +2 $s | cut -f3 | sort -n | head -1); hi=$(tail -n +2 $s | cut -f3 | sort -n | tail -1)
  t0=$(date +%s.%N); $B ld --bfile $R --chr $c --start $lo --end $hi --maf 0 --extract w.ids --out rs 2> rs.err; t1=$(date +%s.%N)
  $P1 --bfile $R --extract w.ids --keep-allele-order --r square bin --threads 48 --out p1 > /dev/null 2>&1; t2=$(date +%s.%N)
  $P2 --bfile $R --extract w.ids --r-unphased square bin --threads 48 --out p2 > /dev/null 2>&1; t3=$(date +%s.%N)
  d=$($PY -I -c "
import numpy as np
p=$p; a=np.fromfile('rs.ld.bin'); b=np.fromfile('p1.ld.bin'); c=np.fromfile('p2.unphased.vcor1.bin')
ok=a.size==p*p==b.size==c.size
print(f'{np.nanmax(np.abs(a-b)):.2e}\t{np.nanmax(np.abs(np.abs(a)-np.abs(c))):.2e}' if ok else 'size\tsize')")
  echo -e "$id\t$p\t$(echo "$t1-$t0"|bc)\t$(echo "$t2-$t1"|bc)\t$(echo "$t3-$t2"|bc)\t$d" >> ../ld_times.tsv
  rm -f rs.ld.bin p1.ld.bin p2.unphased.vcor1.bin
done
touch ../ld.done
