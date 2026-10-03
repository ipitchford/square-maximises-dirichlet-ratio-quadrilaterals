#!/bin/bash
# slot D: wait for the FF replay to finish, then run PP shard 2/3 and SP_col
cd /Users/admin/claude-sandbox/quad-ppw-2026-09-30/darts
until grep -q "^exit" replays/ff.0.log 2>/dev/null; do sleep 15; done
./prod_fam.sh PP 2 3
./prod_fam.sh SP_col 0 1
