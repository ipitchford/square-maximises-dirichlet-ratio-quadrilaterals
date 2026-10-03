#!/bin/bash
# replay_archive.sh STEP...   Replay of this research archive from the directory that contains this script.
#   check      (minutes) manifests; local-step records (verify/xlocal.py + 5 negative controls); dart receipts (verify/xreceipt.py,
#              proof and structure-only, with the structure-only run's runtime anchors); receipt mutation test (verify/test_xreceipt.py);
#              band fail-closed test (corrupted logs must be rejected); python -O run of the record checker
#   local      (about 10 min, 4 GB)      Theorem 2 from the exact assembly          -> local/replay_v2/ (overwritten in place)
#   band       (about 10 min, 3 threads) Lemma B'                                   -> $OUT/band_runs/
#   shell      (about 50 min, 3 threads) Theorem 3, exhaustive                      -> $OUT/shell.out
#   transport  (about 2 h, 3 threads)    the transport cover, exhaustive            -> global2/runs/transport_replay_<date>.jsonl
#   darts      (hours; SHARDS workers, about 1.7 GB each) Theorem 5, proof mode     -> darts/replays/ (new run directory and receipt)
#   all        check local band shell transport darts
# Environment: PY (Python 3.13 with python-flint 0.9.0, numpy, scipy, scikit-fem 12.0.2, mpmath, sympy), JULIA (Julia 1.12.7 with
# Arblib.jl 1.8.1; set JULIA_PROJECT=<this dir>/environment/julia after instantiating it), T (Julia threads, default 3),
# SHARDS (dart workers, default 2), OUT (default ./replay_out).  Exit status 0 iff every requested step passed.
set -uo pipefail
A=$(cd "$(dirname "$0")" && pwd -P)
PY=${PY:-python3}; JULIA=${JULIA:-julia}; T=${T:-3}; SHARDS=${SHARDS:-2}; OUT=${OUT:-$A/replay_out}
mkdir -p "$OUT"
REGION=0,0.225,0,0.225,0,0.225,-0.225,0.225
# The dart receipts shipped with the archive (fresh-extraction replay of 2026-10-02/03; see REPLAY_RECEIPT.md)
PROOF_RECEIPT="replays/receipts/20261002T221832_7c5eb4b8.jsonl"; QUICK_RECEIPT="replays/receipts/20261002T202211_d322f388.jsonl"; PIN_SHA="58bd96bf9ee620e359dd4cdbb420f6c79096bf1cdb1a73e43336647de697659d"
fail=0
ok()  { echo "PASS  $1"; }
bad() { echo "FAIL  $1"; fail=1; }
recorded_root() { "$PY" -c 'import json,sys,os; r=json.loads(open(sys.argv[1]).readline()); print(os.path.dirname(r["data"]["argv"][0]))' "$1"; }
recorded_home() { "$PY" -c 'import json,sys; r=json.loads(open(sys.argv[1]).readline()); print(r["data"]["child_env"]["HOME"])' "$1"; }
anchor() { "$PY" "$A/verify/xreceipt.py" "$1" --pin-sha256 "$PIN_SHA" --mode quick --darts "$A/darts" --recorded-darts "$(recorded_root "$1")" --home "$(recorded_home "$1")" \
           | "$PY" -c 'import json,sys; l=[x for x in sys.stdin if x.startswith("XRECEIPT {")][0]; d=json.loads(l[9:]); print(d["driver_inventory"], d["shard_inventory"])'; }
step_check() {
  (cd "$A" && "$PY" - <<'PYEOF'
import hashlib, sys
bad = 0
for line in open("MANIFEST.sha256"):
    h, f = line.rstrip("\n").split("  ", 1)
    if hashlib.sha256(open(f, "rb").read()).hexdigest() != h: print("MISMATCH", f); bad += 1
sys.exit(1 if bad else 0)
PYEOF
  ) && ok "archive MANIFEST.sha256" || bad "archive MANIFEST.sha256"
  for c in local global2 bandB; do (cd "$A/$c" && "$PY" -c '
import hashlib,sys
b=[f for h,f in (l.rstrip("\n").split("  ",1) for l in open("MANIFEST.sha256") if l.strip()) if hashlib.sha256(open(f,"rb").read()).hexdigest()!=h]
print("\n".join(b)); sys.exit(1 if b else 0)') && ok "$c/MANIFEST.sha256" || bad "$c/MANIFEST.sha256"; done
  "$PY" "$A/verify/xlocal.py" "$A/local/replay_v2/schur_cover2.json" > "$OUT/xlocal.log" 2>&1 && ok "xlocal (local-step records)" || bad "xlocal"
  "$PY" -O "$A/verify/xlocal.py" "$A/local/replay_v2/schur_cover2.json" > "$OUT/xlocal_O.log" 2>&1 && ok "xlocal under python -O" || bad "xlocal -O"
  "$PY" "$A/verify/xlocal_negative_controls.py" "$A/local/replay_v2/schur_cover2.json" > "$OUT/xlocal_negative.log" 2>&1 && ok "xlocal negative controls (5 corruptions rejected)" || bad "xlocal negative controls"
  read DRV SHD < <(anchor "$A/darts/$QUICK_RECEIPT")
  [ -n "${DRV:-}" ] && ok "structure-only dart receipt verified (anchors $DRV $SHD)" || bad "structure-only dart receipt"
  "$PY" "$A/verify/xreceipt.py" "$A/darts/$PROOF_RECEIPT" --pin-sha256 "$PIN_SHA" --mode proof --driver-inv "$DRV" --shard-inv "$SHD" \
        --darts "$A/darts" --recorded-darts "$(recorded_root "$A/darts/$PROOF_RECEIPT")" --home "$(recorded_home "$A/darts/$PROOF_RECEIPT")" > "$OUT/xreceipt_proof.log" 2>&1 \
    && ok "dart proof receipt verified (exact state machine, pinned files, logs, shard outputs)" || bad "dart proof receipt"
  PY="$PY" "$PY" "$A/verify/test_xreceipt.py" "$A/darts/$PROOF_RECEIPT" --pin-sha256 "$PIN_SHA" --mode proof --darts "$A/darts" \
        --recorded-darts "$(recorded_root "$A/darts/$PROOF_RECEIPT")" --home "$(recorded_home "$A/darts/$PROOF_RECEIPT")" > "$OUT/test_xreceipt.log" 2>&1 \
    && ok "receipt mutation test (every semantic mutation rejected)" || bad "receipt mutation test"
  band_failclosed > "$OUT/band_failclosed.log" 2>&1 && ok "band fail-closed test (22 corrupted logs rejected)" || bad "band fail-closed test"
}
band_failclosed() (  # subshell: the cases and expected exits of bandB/replay_failclosed_test.sh, with $PY and $JULIA
  TMP=$(mktemp -d "${TMPDIR:-/tmp}/bandB_failclosed.XXXXXX") || exit 2; nbad=0
  cd "$A/bandB" || exit 2
  for c in ok unknown_status malformed extra_field dup_line dup_path missing_leaf missing_marker fail_floor wrong_box noncanonical bad_rule \
           false_exclusion summary_mismatch summary_extra_field summary_dup_key summary_bad_subtrees summary_bad_depth missing_summary \
           mixed_chart trailing_garbage empty forged_parent_cert; do
    "$PY" corrupt_log.py runs/recert_chart2.jsonl "$TMP/$c.jsonl" "$c" || { nbad=$((nbad+1)); continue; }
    if [ $c = forged_parent_cert ]; then "$JULIA" -t "$T" replay.jl --out "$TMP" --partial "$TMP/$c.jsonl" > "$TMP/$c.out" 2>&1
    else "$JULIA" -t "$T" replay.jl --out "$TMP" --partial --cert-sample 3 "$TMP/$c.jsonl" > "$TMP/$c.out" 2>&1; fi
    got=$?; want=1; [ $c = ok ] && want=3
    [ "$got" = "$want" ] && echo "ok    $c: exit $got" || { echo "WRONG $c: exit $got, expected $want"; nbad=$((nbad+1)); }
  done
  echo "band fail-closed: $([ $nbad = 0 ] && echo PASS || echo "FAIL ($nbad wrong)")"
  case "$TMP" in */bandB_failclosed.??????) rm -rf -- "$TMP" ;; *) echo "refusing to remove $TMP" ;; esac; exit $nbad
)
step_local() { (cd "$A/local" && PY="$PY" bash replay_local_v2.sh > "$OUT/local.log" 2>&1) && ok "local replay (Theorem 2)" || bad "local replay"; }
step_band() { mkdir -p "$OUT/band_runs"; (cd "$A/bandB" && "$JULIA" -t "$T" replay.jl runs/recert_chart1.jsonl runs/recert_chart2.jsonl --out "$OUT/band_runs" > "$OUT/band.log" 2>&1) \
              && grep -q "REPLAY PASSED" "$OUT/band.log" && ok "band replay (Lemma B')" || bad "band replay"; }
step_shell() { (cd "$A/global2" && "$JULIA" -t "$T" replay.jl runs/shell0225_cert.jsonl --region $REGION > "$OUT/shell.log" 2>&1) \
               && grep -q "^REPLAY OK" "$OUT/shell.log" && ok "shell replay (Theorem 3)" || bad "shell replay"; }
step_transport() { R="runs/transport_replay_$(date +%Y%m%dT%H%M%S).jsonl"
  (cd "$A/global2" && "$JULIA" --heap-size-hint=1500M -t "$T" replay_t.jl runs/transport_cert.jsonl --results "$R" > "$OUT/transport.log" 2>&1) \
  && grep -q "REPLAY OK" "$OUT/transport.log" && ok "transport replay ($R)" || bad "transport replay"; }
step_darts() { (cd "$A/darts" && SHARDS="$SHARDS" PY="$(command -v "$PY")" bash replay_all.sh > "$OUT/darts.log" 2>&1) && ok "dart proof replay (Theorem 5)" || bad "dart proof replay"; }
[ $# -eq 0 ] && { sed -n '2,16p' "$0"; exit 2; }
for s in "$@"; do
  case $s in
    check|local|band|shell|transport|darts) "step_$s" ;;
    all) step_check; step_local; step_band; step_shell; step_transport; step_darts ;;
    *) echo "unknown step $s"; exit 2 ;;
  esac
done
exit $fail
