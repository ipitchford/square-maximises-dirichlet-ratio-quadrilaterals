#!/bin/bash
# Re-run of the sampled recompute after the replay3.py rc-shadowing fix (pin 7e564f3f...): every shard below contains FEM leaves.
#   stream A: PP 23/64 (234 FEM of 3,177);  stream B: SP_spk 23/64 (106 FEM), FP 23/64 (16 FEM), SP_col 42/64 (174 FEM of 525).
D=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
PY=/Users/admin/.venvs/analysis-py313/bin/python
OUT=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/final_run/sample2
/bin/mkdir -p "$OUT" || exit 2
PYC=$(/usr/bin/mktemp -d "$OUT/pyc.XXXXXX") || exit 2
one() {   # family tag shard
  local F="$1" L="$2" S="$3"
  /usr/bin/env -i PATH=/usr/bin:/bin HOME=/Users/admin TMPDIR="$OUT" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    DARTS_PIN="$D/review/PINNED.sha256" \
    "$PY" -I -S -B -X "pycache_prefix=$PYC" "$D/bootstrap.py" "$D/replay3.py" "$D/certs4/$L.merged.jsonl" "$F" "$D/manifests/$L.json" \
    --shard "$S/64" --out "$OUT/$L.$S.jsonl" > "$OUT/$L.$S.log" 2>&1
  local rc=$?
  echo "$(/bin/date '+%T') $F $S/64 exit $rc $(/usr/bin/tail -1 "$OUT/$L.$S.log" | /usr/bin/cut -c1-300)"
  return $rc
}
( one PP pp 23 ) > "$OUT/streamA.out" 2>&1 &
A=$!
( rc=0; for X in SP_spk:sp_spk:23 FP:fp:23 SP_col:sp_col:42; do IFS=: read F L S <<< "$X"; one "$F" "$L" "$S" || rc=1; done; exit $rc ) > "$OUT/streamB.out" 2>&1
B=$?
wait $A; A=$?
/bin/cat "$OUT/streamA.out" "$OUT/streamB.out"
echo "SAMPLE2 streamA exit $A streamB exit $B"
[ $A -eq 0 ] && [ $B -eq 0 ]
