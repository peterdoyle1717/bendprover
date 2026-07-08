#!/bin/bash
# Embeddedness pass over the prove run's records: for each completed
# solver chunk NNN.bends, write emb/NNN.tsv of per-net verdicts
# (PASS / PANCAKE / FAIL). Resumable; safe to run while the solver
# fleet is still going (only touches chunks with .done markers).
# No set -e: idle iterations (nothing to do yet) are normal.
#
#   run/doob_emb.sh launch | status
RUN=$HOME/prove_final_v4_50
REPO=$HOME/bendprover
WORKERS=8

loop() {
  mkdir -p $RUN/emb
  while :; do
    todo=$(cd $RUN/out && ls *.done 2>/dev/null | sed 's/.done//' \
           | while read id; do [ -f $RUN/emb/$id.tsv ] || echo $id; done)
    if [ -n "$todo" ]; then
      echo "$todo" | xargs -P $WORKERS -I{} sh -c \
        "nice -19 $REPO/csrc/embcheck_mp $RUN/out/{}.bends > $RUN/emb/{}.tsv.part 2>/dev/null \
         && mv $RUN/emb/{}.tsv.part $RUN/emb/{}.tsv" || true
    fi
    n=$(ls $RUN/emb/*.tsv 2>/dev/null | wc -l)
    echo "$(date +%H:%M) emb chunks $n/960"
    [ "$n" -ge 960 ] && break
    sleep 60
  done
  echo "emb pass complete: $n chunks"
}

launch() {
  nohup bash $REPO/run/doob_emb.sh loop > $RUN/log/emb.log 2>&1 &
  echo "launched pid $!"
}

status() {
  n=$(ls $RUN/emb/*.tsv 2>/dev/null | wc -l)
  awk -F'\t' 'FNR>1 {c[$3]++} END {for (k in c) printf "%s %d\n", k, c[k]}' \
      $RUN/emb/*.tsv 2>/dev/null | sort
  echo "chunks checked $n/960"
  grep -hP "\tFAIL\t" $RUN/emb/*.tsv 2>/dev/null | head -5
}

"$@"
