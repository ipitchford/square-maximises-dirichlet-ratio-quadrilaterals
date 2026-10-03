#!/bin/bash
# Sampled certificate recompute of the re-certified dart logs BEFORE the 4-hour proof replay (coordinator runbook step).
# Runs the replay child exactly as replay_all.py does (venv python -I -S -B, empty pycache prefix, bootstrap.py, sanitised
# environment, DARTS_PIN = the reviewed pin so chart caches are verified before unpickling) on shard 23/64 of every family's merged
# log (the merged logs are written by the QUICK run).  Stream A: PP; stream B: the seven other families in turn.  Exit 0 iff every
# sampled shard exits 0 with failures 0.  Launch under research-run --cpu 2 --memory-gb 3.
D=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
PY=/Users/admin/.venvs/analysis-py313/bin/python
OUT=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/final_run/sample
/bin/mkdir -p "$OUT" || exit 2
PYC=$(/usr/bin/mktemp -d "$OUT/pyc.XXXXXX") || exit 2
one() {   # family tag
  local F="$1" L="$2"
  /usr/bin/env -i PATH=/usr/bin:/bin HOME=/Users/admin TMPDIR="$OUT" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    DARTS_PIN="$D/review/PINNED.sha256" \
    "$PY" -I -S -B -X "pycache_prefix=$PYC" "$D/bootstrap.py" "$D/replay3.py" "$D/certs4/$L.merged.jsonl" "$F" "$D/manifests/$L.json" \
    --shard 23/64 --out "$OUT/$L.23.jsonl" > "$OUT/$L.23.log" 2>&1
  local rc=$?
  echo "$(/bin/date '+%T') $F exit $rc $(/usr/bin/tail -1 "$OUT/$L.23.log" | /usr/bin/cut -c1-260)"
  return $rc
}
( one PP pp ) > "$OUT/streamA.out" 2>&1 &
A=$!
( rc=0; for FL in SP_col:sp_col SP_spk:sp_spk FP:fp PF:pf FF:ff SF_col:sf_col SF_spk:sf_spk; do one "${FL%%:*}" "${FL##*:}" || rc=1; done; exit $rc ) > "$OUT/streamB.out" 2>&1
B=$?
wait $A; A=$?
/bin/cat "$OUT/streamA.out" "$OUT/streamB.out"
echo "SAMPLE streamA exit $A streamB exit $B"
[ $A -eq 0 ] && [ $B -eq 0 ]
