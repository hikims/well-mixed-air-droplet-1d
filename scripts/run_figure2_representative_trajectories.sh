#!/usr/bin/env bash
set -e

# Run representative trajectories for Figure 2.
# Defaults:
#   stationary  s0=0.080
#   oscillatory s0=0.125
#   rotating    s0=0.180

make figure2 USE_TBB="${USE_TBB:-1}" THREADS="${THREADS:-8}" FIG2_TEND="${FIG2_TEND:-120.0}" FIG2_SAVE_EVERY="${FIG2_SAVE_EVERY:-20}"
