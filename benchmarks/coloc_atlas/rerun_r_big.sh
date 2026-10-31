#!/bin/bash
# R rerun for the 3 large tissues: fread decompresses .gz into TMPDIR, and / had only 5.6 GB free, which truncated
# the input in the first run. TMPDIR now on /mnt/hdd2 (NVMe, same class as /).
cd /mnt/hdd2/loc-tran/susie_rs/atlas
export TMPDIR=/mnt/hdd2/loc-tran/susie_rs/atlas/rtmp; mkdir -p $TMPDIR
E=eqtl; PY=/mnt/hdd2/loc-tran/miniconda3/bin/python
cat eqtl/QTD000121* eqtl/QTD000356* eqtl/QTD000116* > /dev/null
: > logs/r_big.log
for t in adipose_visceral:QTD000121 blood:QTD000356 adipose_subcut:QTD000116; do n=${t%%:*}; q=${t##*:}
  /usr/bin/time -f "R_$n wall %e s maxrss %M KB" ../renv/bin/Rscript atlas_r.R loci.tsv $n $E/$q.lbf_variable.txt.gz $E/$q.credible_sets.tsv.gz full/$n.R.tsv 48 2>&1 | grep -E "R atlas|wall|rror" >> logs/r_big.log
  echo "== $n" >> logs/parity_big.log; $PY -I parity.py full/$n.coloc.tsv full/$n.R.tsv >> logs/parity_big.log 2>&1
done
rm -rf $TMPDIR; touch rbig.done
