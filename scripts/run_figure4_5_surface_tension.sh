#!/usr/bin/env bash
set -e

make figure4 USE_TBB="${USE_TBB:-1}" THREADS="${THREADS:-8}" FIG4_TEND="${FIG4_TEND:-200.0}" FIG4_SAVE_EVERY="${FIG4_SAVE_EVERY:-20}"
