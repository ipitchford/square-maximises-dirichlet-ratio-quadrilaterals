#!/bin/bash
# replay_failclosed_test.sh — the replay must FAIL CLOSED: every corrupted copy of a valid log must give exit 1;
# the intact log in test mode gives exit 3 (no failure, not a proof).  Author: Anonymous.
# usage: ./replay_failclosed_test.sh [log]   (default runs/cover_chart2.jsonl)  -> replay_failclosed_test.log
cd "$(dirname "$0")"
SRC=${1:-runs/recert_chart2.jsonl}
TMP=$(mktemp -d "${TMPDIR:-/tmp}/bandB_failclosed.XXXXXX")
[ -n "$TMP" ] && [ -d "$TMP" ] || { echo "no temp dir"; exit 2; }
PY=~/.venvs/analysis-py313/bin/python
JL="nice -n 10 /opt/homebrew/bin/julia -t 3 replay.jl --out $TMP"
nbad=0
run() {  # case expected_exit extra_args...
  local c=$1 want=$2; shift 2
  $PY corrupt_log.py "$SRC" "$TMP/$c.jsonl" "$c" || { echo "corrupt_log failed for $c"; nbad=$((nbad+1)); return; }
  $JL "$@" "$TMP/$c.jsonl" > "$TMP/$c.out" 2>&1; local got=$?
  local why=$(grep -m1 "^FAIL" "$TMP/$c.out" | cut -c1-150)
  if [ "$got" = "$want" ]; then echo "ok    $c: exit $got   $why"; else echo "WRONG $c: exit $got, expected $want   $why"; nbad=$((nbad+1)); fi
}
echo "fail-closed tests on $SRC  $(date)"
run ok 3 --partial --cert-sample 3
for c in unknown_status malformed extra_field dup_line dup_path missing_leaf missing_marker fail_floor wrong_box \
         noncanonical bad_rule false_exclusion summary_mismatch summary_extra_field summary_dup_key summary_bad_subtrees \
         summary_bad_depth missing_summary mixed_chart trailing_garbage empty; do
  run $c 1 --partial --cert-sample 3
done
run forged_parent_cert 1 --partial          # all certificates re-run: only the parent's certificate can fail here
run ok 1 --cert-sample 3                    # without --partial both charts are required
echo "fail-closed test: $([ $nbad = 0 ] && echo PASS || echo "FAIL ($nbad wrong)")  $(date)"
case "$TMP" in */bandB_failclosed.??????) rm -rf -- "$TMP" ;; *) echo "refusing to remove $TMP" ;; esac
exit $nbad
