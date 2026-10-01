#!/usr/bin/env bash
# Runs a command and samples its CPU time every INTERVAL seconds (default 0.25).
# Reports average and peak cores busy, and peak thread count.
# Usage: cpu_monitor.sh [-i INTERVAL] -- command args...
# Works on Linux (/proc) and macOS (ps).
set -euo pipefail
interval=0.25
if [ "${1:-}" = "-i" ]; then interval=$2; shift 2; fi
[ "${1:-}" = "--" ] && shift

"$@" &
pid=$!

cpu_seconds() {  # cumulative user+sys CPU seconds of $pid
  if [ -r "/proc/$pid/stat" ]; then
    awk -v hz="$(getconf CLK_TCK)" '{print ($14 + $15) / hz}' "/proc/$pid/stat" 2>/dev/null
  else
    ps -o time= -p "$pid" 2>/dev/null | awk -F: '{ if (NF == 3) s = $1*3600 + $2*60 + $3; else s = $1*60 + $2; print s }'
  fi
}
now_s() {  # wall clock in seconds with sub-second resolution
  if [ -r /proc/uptime ]; then awk '{print $1}' /proc/uptime
  else perl -MTime::HiRes=time -e 'printf "%.3f", time'; fi
}
threads() {
  if [ -r "/proc/$pid/status" ]; then awk '/^Threads:/{print $2}' "/proc/$pid/status" 2>/dev/null
  else ps -M -p "$pid" 2>/dev/null | tail -n +2 | wc -l | tr -d ' '; fi
}

t0=$(now_s)
prev_t=$t0; prev_c=0; peak=0; peak_at=0; max_thr=0; samples=""
while kill -0 "$pid" 2>/dev/null; do
  sleep "$interval"
  c=$(cpu_seconds) || break; [ -z "$c" ] && break
  n=$(threads); [ -n "$n" ] && [ "$n" -gt "$max_thr" ] && max_thr=$n
  now=$(now_s)
  cores=$(awk -v c="$c" -v pc="$prev_c" -v t="$now" -v pt="$prev_t" 'BEGIN{ d = t - pt; printf "%.2f", (d > 0 ? (c - pc) / d : 0) }')
  el=$(awk -v t="$now" -v t0="$t0" 'BEGIN{printf "%.2f", t - t0}')
  samples="$samples $el:$cores"
  if awk -v a="$cores" -v b="$peak" 'BEGIN{exit !(a > b)}'; then peak=$cores; peak_at=$el; fi
  prev_t=$now; prev_c=$c
done
wait "$pid" || true
end=$(now_s)
total=$(awk -v c="$prev_c" -v e="$end" -v t0="$t0" 'BEGIN{printf "%.2f %.2f %.2f", e - t0, c, c / (e - t0)}')
set -- $total
echo "wall ${1}s  cpu ${2}s  avg ${3} cores  peak ${peak} cores (at ${peak_at}s)  max threads ${max_thr}" >&2
[ -n "${CPU_MONITOR_TRACE:-}" ] && echo "trace:$samples" >&2
exit 0
