#!/bin/bash
# Queued launch of the proof replay through research-run, after the re-certification is complete and the round-4 review has
# cleared.  Started by the coordinator (not by the agent that prepared it), never nested inside research-run.
# Sizing rule (coordinator, 2026-10-02): SHARDS = largest <= 4 with SHARDS x (1.50 GB x L / 621,307) x 1.1 + 0.4 <= 6.0 GB, where
# 1.50 GB is the measured peak RSS (/usr/bin/time -l) of one replay3 shard of certs3/pp.merged.jsonl (621,307 lines; every shard
# parses the whole log) and L is the line count of the certs4 PP worker logs (= the merged log the PP shards parse); the driver
# replays the logs one after another, so the PP shards are the peak.  --memory-gb = that figure rounded up to 0.5 GB.
# The driver aborts within seconds (certs_complete) if certs4 is incomplete or has a FAIL leaf.
D=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
set -o pipefail
L=$(/bin/cat "$D"/certs4/pp.w*.jsonl | /usr/bin/wc -l | /usr/bin/tr -d ' ') || { echo "no certs4 PP worker logs"; exit 2; }
read -r SH MEM < <(/usr/bin/awk -v L="$L" 'BEGIN { per = 1.50 * L / 621307 * 1.1; s = 0
  for (k = 4; k >= 1; k--) if (k * per + 0.4 <= 6.0) { s = k; break }
  m = s * per + 0.4; h = int(m * 2); if (h < m * 2) h++; printf "%d %.1f\n", s, h / 2 }')
if ! [[ "$SH" =~ ^[1-4]$ ]]; then echo "launch_replay: PP log of $L lines does not fit the 6 GB rule"; exit 2; fi
{
  echo "launch_replay: $(/bin/date '+%F %T') certs4 PP lines $L -> SHARDS=$SH memory-gb=$MEM"
  /Users/admin/.local/bin/research-run run --cpu "$SH" --memory-gb "$MEM" --max-wait 0 --label darts-replay -- \
    /usr/bin/env PIN="$D/review/PINNED.sha256" SHARDS="$SH" /bin/bash "$D/replay_all.sh"
  echo "launch_replay: research-run exit $? $(/bin/date '+%F %T')"
} >> "$D/certs4/launch.out" 2>&1
