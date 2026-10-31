#!/bin/bash
# Full atlas, 6 GTEx tissues. Files are read once into the page cache before timing so both tools start warm.
cd /mnt/hdd2/loc-tran/susie_rs/atlas
until [ -f eqtl/fetch.done ]; do sleep 30; done
B=../susie-rs/target/release/susiers; E=eqtl; PY=/mnt/hdd2/loc-tran/miniconda3/bin/python
T="uterus:QTD000346 ovary:QTD000291 vagina:QTD000351 adipose_visceral:QTD000121 blood:QTD000356 adipose_subcut:QTD000116"
Q=""; for t in $T; do n=${t%%:*}; q=${t##*:}; Q="$Q --qtl $n=$E/$q.lbf_variable.txt.gz,$E/$q.credible_sets.tsv.gz"; done
cat eqtl/*.gz > /dev/null
mkdir -p full logs
/usr/bin/time -f "rust_all wall %e s maxrss %M KB" $B atlas --loci loci.tsv $Q --out full > logs/rust.log 2>&1
for t in $T; do n=${t%%:*}; q=${t##*:}
  /usr/bin/time -f "R_$n wall %e s maxrss %M KB" ../renv/bin/Rscript atlas_r.R loci.tsv $n $E/$q.lbf_variable.txt.gz $E/$q.credible_sets.tsv.gz full/$n.R.tsv 48 2>&1 | grep -E "R atlas|wall" >> logs/r.log
done
/usr/bin/time -f "rust_all_rep wall %e s maxrss %M KB" $B atlas --loci loci.tsv $Q --out full_rep > logs/rust_rep.log 2>&1
for t in $T; do n=${t%%:*}; echo "== $n"; $PY -I parity.py full/$n.coloc.tsv full/$n.R.tsv; done > logs/parity.log 2>&1
touch atlas.done
