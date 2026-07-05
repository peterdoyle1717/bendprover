#!/bin/bash
# doob production run: all prime nets v4..50, seeded from the May-2026
# pipeline's LM obj solutions. Stage once, then chunked workers
# (1000 nets per python process, nice 19, $WORKERS cores).
#
#   ./run/doob_v4_50.sh stage    # names, worklist, chunks, alignment checks
#   ./run/doob_v4_50.sh launch   # start workers under nohup
#   ./run/doob_v4_50.sh status   # chunk + pass/fail tally
set -e
PIPE=$HOME/neo/data/pipeline_v4_50_20260514T103000Z
CLERS_BIN=$HOME/neo/clers/bin/clers
RUN=$HOME/bendprover_run_v4_50
REPO=$HOME/bendprover
WORKERS=90

stage() {
  # identity source: solvermanifest/all.tsv = index, status, netcode, objpath
  # (the solver wrote each obj FROM that netcode; netcodes.txt/objfiles.txt
  # are NOT line-aligned with each other)
  MAN=$PIPE/euclid/solvermanifest/all.tsv
  mkdir -p $RUN/chunks $RUN/out $RUN/log
  cd $RUN
  wc -l $MAN
  echo "status tally:"; cut -f2 $MAN | sort | uniq -c
  echo "naming netcodes (one clers process)..."
  cut -f3 $MAN | $CLERS_BIN name > names.txt
  cut -f4 $MAN > objs.txt
  wc -l names.txt objs.txt
  echo "unique names: $(sort -u names.txt | wc -l)"
  paste names.txt objs.txt > worklist.tsv
  # alignment spot-check: obj faces must equal the manifest netcode
  for n in 1 173000 4200000 8239684; do
    row=$(awk -v n=$n 'NR==n' $MAN)
    nc=$(printf "%s" "$row" | cut -f3); obj=$(printf "%s" "$row" | cut -f4)
    of=$(grep '^f ' $obj | awk '{printf "%s%d,%d,%d",sep,$2,$3,$4; sep=";"} END{print ""}')
    [ "$nc" = "$of" ] && echo "align ok line $n" || { echo "ALIGN FAIL line $n"; exit 1; }
  done
  split -d -a 5 -l 1000 worklist.tsv chunks/chunk_
  ls chunks | wc -l
}

worker() {  # $1 = chunk file
  id=$(basename $1 | sed s/chunk_//)
  [ -f $RUN/out/$id.done ] && exit 0
  CLERS_BIN=$CLERS_BIN python3 $REPO/bendprover.py --batch $1 $RUN/out/$id \
      > $RUN/log/$id.log 2>&1 && touch $RUN/out/$id.done
}

launch() {
  cd $RUN
  ls chunks/chunk_* | nohup nice -19 xargs -P $WORKERS -I{} bash $REPO/run/doob_v4_50.sh worker {} \
      > log/xargs.log 2>&1 &
  echo "launched pid $!"
}

status() {
  total=$(ls $RUN/chunks | wc -l)
  done_=$(ls $RUN/out/*.done 2>/dev/null | wc -l)
  pass=$(cat $RUN/log/*.log 2>/dev/null | grep -c '^PASS' || true)
  fail=$(cat $RUN/log/*.log 2>/dev/null | grep -c '^FAIL' || true)
  echo "chunks $done_/$total   nets: $pass pass, $fail fail"
}

"$@"
