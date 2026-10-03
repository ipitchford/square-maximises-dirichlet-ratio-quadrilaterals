#!/bin/bash
# Keepalive for pool workers 1-6 (worker 7 belongs to pool_scaler.sh): every 5 min, a worker whose lock pid is dead and whose
# output does not end with DONE is restarted (it resumes its claimed items from its own log).  Log: pool_runs/keepalive.log.
cd "$(dirname "$0")"
LOG=pool_runs/keepalive.log
while true; do
  for i in 1 2 3 4 5 6; do
    p=$(cat pool_runs/w$i.lock 2>/dev/null)
    if [ -n "$p" ] && kill -0 "$p" 2>/dev/null; then continue; fi
    if tail -1 pool_runs/w$i.out 2>/dev/null | grep -q '^DONE'; then continue; fi
    echo "=== restart $(date +%H:%M:%S) (keepalive: worker $i was not running)" >> pool_runs/w$i.err
    ./pool_launch.sh $i $i >> $LOG 2>&1
    echo "$(date '+%F %T') restarted worker $i" >> $LOG
  done
  sleep 300
done
