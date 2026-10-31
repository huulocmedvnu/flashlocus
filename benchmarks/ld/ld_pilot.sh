#!/bin/bash
# LD pilot: same SNP list per window, 48 threads for every tool, wall time per invocation (bim parsing included for all).
cd /mnt/hdd2/loc-tran/susie_rs; mkdir -p ldbench/pilot; cd ldbench/pilot
R=/mnt/hdd2/loc-tran/infertility_work/ref/g1000_eur
B=../../susie-rs/target/release/susiers; P1=/mnt/hdd2/loc-tran/icp_mr/tools/plink; P2=../../tools/plink2/plink2
cat $R.bed $R.bim > /dev/null
for id in $(awk -F'\t' 'NR>1{print $1"\t"$8}' ../../real/*/prep_summary.tsv | sort -k2,2nr | head -3 | cut -f1); do
  t=${id//[0-9]/}; s=../../real/$t/$id.snps.tsv
  tail -n +2 $s | cut -f1 > $id.ids; c=$(sed -n 2p $s | cut -f2); lo=$(tail -n +2 $s | cut -f3 | sort -n | head -1); hi=$(tail -n +2 $s | cut -f3 | sort -n | tail -1)
  echo "== $id p=$(wc -l < $id.ids)"
  /usr/bin/time -f "susiers %e s" $B ld --bfile $R --chr $c --start $lo --end $hi --maf 0 --extract $id.ids --out $id.rs 2>&1 | grep -E "s$"
  /usr/bin/time -f "plink1.9 %e s" $P1 --bfile $R --extract $id.ids --keep-allele-order --r square bin --threads 48 --out $id.p1 > /dev/null 2>&1 ; tail -1 $id.p1.log >/dev/null
  grep -h "plink1.9" /dev/null
  /usr/bin/time -f "plink2 %e s" $P2 --bfile $R --extract $id.ids --r-unphased square bin --threads 48 --out $id.p2 > $id.p2.stdout 2>&1; grep -E "s$" $id.p2.stdout | tail -1
done 2>&1
