#!/bin/bash
# run_cover.sh — the rigorous cover of both charts (BAND_SPEC §4), sequential, at nice 10.  Author: Anonymous.
# Fails (non-zero exit) as soon as any cover run fails.  Threads: 3 (coordinator limit, 2026-10-01 afternoon).
set -euo pipefail
cd "$(dirname "$0")"
T=${T:-3}
for chart in 2 1; do
  echo "chart $chart: $(date) threads $T"
  if ! nice -n 10 /opt/homebrew/bin/julia -t "$T" cover.jl "$chart" "runs/cover_chart$chart.jsonl" > "runs/cover_chart$chart.out" 2>&1; then
    echo "chart $chart FAILED $(date)"; exit 1
  fi
  echo "chart $chart exit 0 $(date)"
done
