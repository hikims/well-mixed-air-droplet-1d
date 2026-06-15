#!/usr/bin/env bash
set -e

# Run a long phase sweep in the background on macOS.
# Usage:
#   bash scripts/run_phase_long_background.sh

mkdir -p logs results figure

LOG="logs/phase_long_$(date +%Y%m%d_%H%M%S).log"

echo "Starting long phase sweep..."
echo "Log file: $LOG"

nohup caffeinate -dimsu make phase_long_all USE_TBB=1 THREADS="${THREADS:-8}" \
    > "$LOG" 2>&1 &

echo $! > logs/phase_long.pid

echo "PID: $(cat logs/phase_long.pid)"
echo "Progress:"
echo "  tail -f $LOG"
