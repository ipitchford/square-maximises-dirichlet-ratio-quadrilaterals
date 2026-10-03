#!/bin/bash
# Queued launch of the re-certification through research-run.  Started by the coordinator (not by the agent that prepared it),
# never nested inside research-run.  research-run execs the job as /usr/bin/nice -n 10 CMD; recert_all.sh keeps the TARGET niceness 10 (no double nice).
# Sizing (round 4, measured with /usr/bin/time -l): one recert worker on heavy PP roots (BD k=6, FARL k=5, INC; 2,148 leaves)
# peaks at 0.43 GB RSS; resume bookkeeping of a PP share adds <= 0.07 GB (status-only, measured with tracemalloc) and the FARL
# problem cache ~0.02 GB: 0.52 GB per worker.  Rule: W = largest <= 4 with W x 0.52 x 1.2 + 0.3 <= 2.5 GB  ->  W = 3;
# --memory-gb = 3 x 0.52 x 1.2 + 0.3 = 2.17 -> 2.5.  The job's exit status is recert_all.sh's: 0 iff certs_complete.py finds
# every root of every family done with no FAIL leaf.
D=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
W=3; MEM=2.5
/bin/mkdir -p "$D/certs4" || exit 2
{
  echo "launch_recert: $(/bin/date '+%F %T') W=$W memory-gb=$MEM"
  /Users/admin/.local/bin/research-run run --cpu "$W" --memory-gb "$MEM" --max-wait 0 --label darts-recert -- /bin/bash "$D/recert_all.sh" "$W"
  echo "launch_recert: research-run exit $? $(/bin/date '+%F %T')"
} >> "$D/certs4/launch.out" 2>&1
