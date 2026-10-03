#!/bin/bash
# Re-certification of all cover5 family trees with the repaired verifier (recert.py): N workers (default 6 = the queue's CPU
# reservation), each working through the families in a fixed order (largest first); a root is claimed by exactly one worker
# (<dir>/claims/<fam>/<root>); rerunning this script with the SAME N resumes (claims and per-worker outputs are kept).
#   usage: recert_all.sh [N]          N canonical integer 1..64 (default 6)
#   output: <dir>/<fam>.w<ID>.jsonl, progress <dir>/recert.w<ID>.out / .err; <dir> = certs4 (RECERT_DIR: test directories only)
# Launch-path hardening (round 4): the package directory and the interpreter are ABSOLUTE constants (never derived from $0 or
# PATH); every tool is called by absolute path; the workers run isolated (-I -S -B, empty pycache prefix, bootstrap.py) under a
# sanitised environment (env -i); the niceness TARGET is 10 (the increment is 10 minus the current niceness, never below 0, so a
# launcher that already applies nice 10, e.g. research-run, gives 10 and not 20).  The first start records N and the code
# identity (recert.py --code-id) in <dir>/RECERT_N and <dir>/RECERT_CODE; a resume with another N or other code is REFUSED
# (a different N would orphan the claims of the missing worker IDs).  When all workers have finished, certs_complete.py checks
# that every root is done (exit status of this script = its verdict; any FAIL leaf is listed).
# Afterwards: PIN=review/PINNED.sha256 SHARDS=6 ./replay_all.sh
D=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
PY=/Users/admin/.venvs/analysis-py313/bin/python
ALLF="PP SP_spk SP_col FP PF FF SF_spk SF_col"
N="${1:-6}"
if ! [[ "$N" =~ ^[1-9][0-9]?$ ]] || [ "$N" -gt 64 ]; then echo "N must be a canonical integer in 1..64, got '$N'"; exit 2; fi
OUT="${RECERT_DIR:-$D/certs4}"
FAMS="$ALLF"
if [ "$OUT" != "$D/certs4" ]; then                                    # test directories only: under the package, absolute
  case "$OUT" in "$D"/recert_tests/*|"$D"/replays/test_runs/*) ;; *) echo "RECERT_DIR must be $D/certs4 or a test directory under recert_tests/ or replays/test_runs/"; exit 2 ;; esac
  if [ -n "${RECERT_FAMILIES-}" ]; then
    for F in $RECERT_FAMILIES; do case " $ALLF " in *" $F "*) ;; *) echo "unknown family '$F'"; exit 2 ;; esac; done
    FAMS="$RECERT_FAMILIES"
  fi
elif [ -n "${RECERT_FAMILIES-}" ]; then echo "RECERT_FAMILIES is allowed for test directories only"; exit 2
fi
CUR=$(/bin/ps -o nice= -p $$ | /usr/bin/tr -d ' ')
if ! [[ "$CUR" =~ ^-?[0-9]+$ ]]; then echo "cannot read the current niceness ('$CUR')"; exit 2; fi
INC=$(( CUR < 10 ? 10 - CUR : 0 ))
/bin/mkdir -p "$OUT" "$D/replays" || exit 2
PYC=$(/usr/bin/mktemp -d "$D/replays/.pycache.XXXXXX") || exit 2
py() { /usr/bin/env -i PATH=/usr/bin:/bin HOME="$HOME" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 RECERT_DIR="$OUT" \
       /usr/bin/nice -n "$INC" "$PY" -I -S -B -X "pycache_prefix=$PYC" "$D/bootstrap.py" "$@"; }
CODE=$(py "$D/recert.py" --code-id) || { echo "recert.py --code-id failed"; exit 2; }
for KV in "RECERT_N=$N" "RECERT_CODE=$CODE" "RECERT_FAMILIES=$FAMS"; do
  F="$OUT/${KV%%=*}"; V="${KV#*=}"
  if [ -e "$F" ]; then
    if [ "$(/bin/cat "$F")" != "$V" ]; then echo "$F records '$(/bin/cat "$F")' but this start has '$V': refusing to resume"; exit 2; fi
  else
    ( set -o noclobber; echo "$V" > "$F" ) || { echo "cannot create $F"; exit 2; }
  fi
done
echo "recert_all: N=$N dir=$OUT families='$FAMS' niceness $CUR -> +$INC code $CODE start $(/bin/date '+%F %T')"
for F in $FAMS; do py "$D/recert.py" --split "$F" > /dev/null || { echo "split $F failed"; exit 1; }; done
for ((w=1; w<=N; w++)); do
  ( for F in $FAMS; do
      py "$D/recert.py" "$F" "$w" >> "$OUT/recert.w$w.out" 2>> "$OUT/recert.w$w.err" || echo "worker $w family $F exit $?" >> "$OUT/recert.w$w.err"
    done; echo "ALLDONE $(/bin/date '+%F %T')" >> "$OUT/recert.w$w.out" ) &
done
wait
echo "recert_all: workers finished $(/bin/date '+%F %T'); completeness check:"
py "$D/certs_complete.py" "$OUT" $FAMS
RC=$?
/bin/rmdir "$PYC" 2>/dev/null
exit $RC
