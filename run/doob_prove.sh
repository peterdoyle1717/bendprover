#!/bin/bash
# The production run: every prime net v4..50 solved and certified from
# scratch by euclid_lm_mp --prove, records in the adopted format
# (FORMAT.md), names from the canonical CLERS lists.
#
#   run/doob_prove.sh stage | launch | status | manifest
set -e
NETCODES=$HOME/neo/data/pipeline_v4_50_20260514T103000Z/netcodes.txt
CLERS=$HOME/neo/clers/bin/clers
RUN=$HOME/prove_final_v4_50
REPO=$HOME/bendprover
WORKERS=96
CHUNKS=960

stage() {
  mkdir -p $RUN/chunks $RUN/out $RUN/log
  echo "naming $(wc -l < $NETCODES) netcodes (one clers process)..."
  $CLERS name < $NETCODES > $RUN/names.txt
  paste -d' ' $RUN/names.txt $NETCODES > $RUN/input.txt
  split -d -a 3 -n l/$CHUNKS $RUN/input.txt $RUN/chunks/chunk_
  ls $RUN/chunks | wc -l
}

worker() {
  id=$(basename $1 | sed s/chunk_//)
  [ -f $RUN/out/$id.done ] && exit 0
  $REPO/csrc/euclid_lm_mp --prove --batch < $1 \
      > $RUN/out/$id.bends 2> $RUN/log/$id.err \
    && touch $RUN/out/$id.done
}

launch() {
  ls $RUN/chunks/chunk_* | nohup nice -19 xargs -P $WORKERS -I{} \
      bash $REPO/run/doob_prove.sh worker {} > $RUN/log/xargs.log 2>&1 &
  echo "launched pid $!"
}

status() {
  done_=$(ls $RUN/out/*.done 2>/dev/null | wc -l)
  recs=$(cat $RUN/out/*.bends 2>/dev/null | grep -c "^end" || true)
  fails=$(cat $RUN/out/*.bends 2>/dev/null | grep -c "^# failed" || true)
  echo "chunks done $done_/$CHUNKS  records $recs  failed $fails"
}

manifest() {
  cd $RUN/out
  for f in *.bends; do
    n=$(grep -c "^end" $f)
    mx=$(grep "^benderr" $f | awk '{if ($2+0 > m) m = $2+0} END {print m+0}')
    echo "file=$f  records=$n  max_benderr=$mx  sha256=$(sha256sum $f | cut -d\" \" -f1)  date=$(date +%F)"
  done > $RUN/MANIFEST.txt
  wc -l $RUN/MANIFEST.txt
}

"$@"
