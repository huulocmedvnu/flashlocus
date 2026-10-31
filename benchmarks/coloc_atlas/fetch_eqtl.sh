#!/bin/bash
# eQTL Catalogue GTEx SuSiE files, each split into 24 byte ranges fetched in parallel (EBI ~0.4 MB/s per connection)
set -e
D=${BENCH_ROOT:-.}/atlas/eqtl; mkdir -p $D; cd $D
base=http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015
get() {  # url out
  u=$1; o=$2; [ -f $o.ok ] && return
  s=$(curl -sI $u | grep -i content-length | tr -d '\r' | awk '{print $2}')
  k=24; c=$(( (s + k - 1) / k ))
  for i in $(seq 0 $((k-1))); do a=$((i*c)); b=$(( (i+1)*c - 1 )); [ $b -ge $s ] && b=$((s-1))
    (until curl -s --retry 5 -r $a-$b -o $o.part$i $u && [ $(stat -c%s $o.part$i) -eq $((b-a+1)) ]; do sleep 5; done) & done
  wait; cat $(for i in $(seq 0 $((k-1))); do echo $o.part$i; done) > $o; rm $o.part*
  [ $(stat -c%s $o) -eq $s ] && gzip -t $o && touch $o.ok
}
for q in QTD000346 QTD000291 QTD000351 QTD000121 QTD000356 QTD000116; do
  get $base/$q/$q.credible_sets.tsv.gz $q.credible_sets.tsv.gz
  get $base/$q/$q.lbf_variable.txt.gz $q.lbf_variable.txt.gz
  echo "$q done $(date)"
done
touch fetch.done
