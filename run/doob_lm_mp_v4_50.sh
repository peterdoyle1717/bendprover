#!/bin/bash
# doob production: all prime nets v4..50 solved from scratch by
# euclid_lm_mp (128-bit gated LM, from the wish), bends-only output at
# the computed precision. 960 chunks for restart granularity, 96
# nice-19 workers, one solver process per chunk.
#
#   run/doob_lm_mp_v4_50.sh stage | launch | status
set -e
NETCODES=$HOME/neo/data/pipeline_v4_50_20260514T103000Z/netcodes.txt
RUN=$HOME/lm_mp_run_v4_50
REPO=$HOME/bendprover
WORKERS=96
CHUNKS=960

stage() {
  mkdir -p $RUN/chunks $RUN/out $RUN/log
  wc -l $NETCODES
  split -d -a 3 -n l/$CHUNKS $NETCODES $RUN/chunks/chunk_
  ls $RUN/chunks | wc -l
}

worker() {
  id=$(basename $1 | sed s/chunk_//)
  [ -f $RUN/out/$id.done ] && exit 0
  $REPO/csrc/euclid_lm_mp --bends-only --batch < $1 \
      > $RUN/out/$id.bends 2> $RUN/log/$id.err \
    && touch $RUN/out/$id.done
}

launch() {
  ls $RUN/chunks/chunk_* | nohup nice -19 xargs -P $WORKERS -I{} \
      bash $REPO/run/doob_lm_mp_v4_50.sh worker {} > $RUN/log/xargs.log 2>&1 &
  echo "launched pid $!"
}

status() {
  done_=$(ls $RUN/out/*.done 2>/dev/null | wc -l)
  ok=$(cat $RUN/out/*.bends 2>/dev/null | grep -c "^=== .* ok" || true)
  fail=$(cat $RUN/out/*.bends 2>/dev/null | grep -c "^=== .* fail" || true)
  echo "chunks done $done_/$CHUNKS   nets ok $ok   fail $fail"
}

"$@"
