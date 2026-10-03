#!/bin/bash
# Pool scaler (coordinator's rule): every 5 min, check that no ffmpeg / chrome-headless process runs, the 5-minute load average is
# below 30 and free + inactive memory is at least 2 GB.  After three consecutive good checks (10 min) start pool worker 7
# (nice 10).  While worker 7 runs, stop it again if free + inactive memory falls below 2 GB (it resumes later from its log).
# Log: pool_runs/scaler.log.  Usage: nohup ./pool_scaler.sh &
cd "$(dirname "$0")"
LOG=pool_runs/scaler.log
good=0
mem_gb() { vm_stat | awk '/page size of/ {ps=$8} /Pages free/ {f=$3} /Pages inactive/ {i=$3} END {gsub("\\.","",f); gsub("\\.","",i); printf "%.2f", (f+i)*ps/1073741824}'; }
w7_pid() { [ -f pool_runs/w7.lock ] && p=$(cat pool_runs/w7.lock) && kill -0 "$p" 2>/dev/null && echo "$p"; }
while true; do
  load5=$(sysctl -n vm.loadavg | awk '{print $3}')
  mem=$(mem_gb)
  busy=$(pgrep -f "ffmpeg|chrome-headless" | wc -l | tr -d ' ')
  pid7=$(w7_pid)
  if [ -n "$pid7" ]; then
    if awk -v m="$mem" 'BEGIN {exit !(m < 2.0)}'; then
      kill "$pid7"; echo "$(date '+%F %T') memory ${mem} GB < 2 GB: stopped worker 7 (pid $pid7)" >> $LOG; good=0
    fi
  else
    if [ "$busy" = "0" ] && awk -v l="$load5" -v m="$mem" 'BEGIN {exit !(l < 30 && m >= 2.0)}'; then good=$((good + 1)); else good=0; fi
    echo "$(date '+%F %T') load5 $load5 mem ${mem} GB render-procs $busy good-checks $good" >> $LOG
    if [ $good -ge 3 ]; then
      echo "=== start $(date +%H:%M:%S) (scaler: 7th worker)" >> pool_runs/w7.err
      ./pool_launch.sh 7 7 >> $LOG 2>&1; echo "$(date '+%F %T') started worker 7" >> $LOG; good=0
    fi
  fi
  sleep 300
done
