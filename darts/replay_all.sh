#!/bin/bash
# Launcher of the aggregate fail-closed proof replay replay_all.py (see its docstring).  The package directory is the directory of
# the physical launcher file, symlinks resolved (portability patch of 2026-10-02; before it, an absolute constant of the producer's
# machine), so a symlink to this launcher still selects the canonical driver; the
# interpreter is PY (an absolute path to an executable; default: the producer's interpreter).  The driver
# runs isolated (-I), without the site module (-S, no *.pth startup hooks), without bytecode (-B and a new empty pycache prefix),
# under a sanitised environment (env -i: only PATH=/usr/bin:/bin, HOME and the BLAS thread variables).  Environment of the caller:
#   SHARDS  canonical integer 1..64 (default 6 = the queue's CPU reservation)        QUICK  unset/empty (proof) or exactly 1 (structure only, never PASS)
#   PIN     the reviewed pin list (default <package>/review/PINNED.sha256)   CERTS_DIR  re-certified family logs (default certs4)
#   PY      absolute path of the Python interpreter (default /Users/admin/.venvs/analysis-py313/bin/python)
# Any other value is rejected (exit 2).  Exit 0 iff the run is a proof PASS.  Extra arguments are passed to the driver (machinery
# tests: --mode test --only TAGS --gates none --run-id ID; a test run can never PASS).
SELF=$(/usr/bin/readlink -f "$0") || exit 2                       # the physical launcher file (a symlink to it resolves here)
D=$(/usr/bin/dirname "$SELF") || exit 2
PY="${PY:-/Users/admin/.venvs/analysis-py313/bin/python}"
case "$PY" in /*) ;; *) echo "PY must be an absolute path, got '$PY'"; exit 2 ;; esac
[ -x "$PY" ] || { echo "PY is not an executable file: $PY"; exit 2; }
SH="${SHARDS:-6}"
if ! [[ "$SH" =~ ^[1-9][0-9]?$ ]] || [ "$SH" -gt 64 ]; then echo "SHARDS must be a canonical integer in 1..64, got '$SH'"; exit 2; fi
case "${QUICK-}" in
  "") MODE=proof ;;
  1) MODE=quick ;;
  *) echo "QUICK must be unset or exactly 1, got '$QUICK'"; exit 2 ;;
esac
PIN_ARG="${PIN:-$D/review/PINNED.sha256}"
CERTS_ARG="${CERTS_DIR:-certs4}"
/bin/mkdir -p "$D/replays" || exit 2
PYC=$(/usr/bin/mktemp -d "$D/replays/.pycache.XXXXXX") || exit 2
exec /usr/bin/env -i PATH=/usr/bin:/bin HOME="$HOME" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  "$PY" -I -S -B -X "pycache_prefix=$PYC" "$D/bootstrap.py" "$D/replay_all.py" --pin "$PIN_ARG" --shards "$SH" --mode "$MODE" --certs "$CERTS_ARG" "$@"
