#!/bin/bash
# Dry runs of the exact launch path under the real launcher, 1 CPU each.  NOT run by the agent that prepared them: executing
# research-run was denied to it by the permission system; whether these are run is the user's decision.
#   smoke_launch.sh recert   re-certify family FF (42 roots, 540 leaves) with recert_all.sh 1 into recert_tests/smoke_<stamp>
#   smoke_launch.sh replay   replay crit_collar (4 legacy certificates) and the newest smoke FF directory through the full
#                            proof-path harness in test mode (pin, completeness, merge, binding, FP checks, receipt; never PASS:
#                            expected exit 3 = "MACHINERY TEST OK (NOT A PROOF)"), SHARDS=1
# Output appended to recert_tests/smoke.out.  Expected wall time: about 1 minute each once admitted.
D=/Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
RR=/Users/admin/.local/bin/research-run
case "$1" in
  recert)
    S="$D/recert_tests/smoke_$(/bin/date '+%Y%m%dT%H%M%S')"
    { echo "smoke recert $(/bin/date '+%F %T') dir $S"
      "$RR" run --cpu 1 --memory-gb 1 --max-wait 0 --label darts-smoke-recert -- \
        /usr/bin/env RECERT_DIR="$S" RECERT_FAMILIES=FF /bin/bash "$D/recert_all.sh" 1
      echo "smoke recert: research-run exit $? $(/bin/date '+%F %T')"; } >> "$D/recert_tests/smoke.out" 2>&1 ;;
  replay)
    S=$(/bin/ls -d "$D"/recert_tests/smoke_2* 2>/dev/null | /usr/bin/tail -1)
    [ -n "$S" ] || { echo "run smoke_launch.sh recert first"; exit 2; }
    { echo "smoke replay $(/bin/date '+%F %T') certs ${S#$D/}"
      "$RR" run --cpu 1 --memory-gb 1 --max-wait 0 --label darts-smoke-replay -- \
        /usr/bin/env PIN="$D/review/PINNED.sha256" SHARDS=1 /bin/bash "$D/replay_all.sh" --mode test --only crit_collar,ff --gates none --certs "${S#$D/}"
      echo "smoke replay: research-run exit $? $(/bin/date '+%F %T')"; } >> "$D/recert_tests/smoke.out" 2>&1 ;;
  *) echo "usage: smoke_launch.sh recert|replay"; exit 2 ;;
esac
