#!/usr/bin/env bash
set -e

make figure6 USE_TBB="${USE_TBB:-1}" THREADS="${THREADS:-8}" \
    FIG6_SMIN="${FIG6_SMIN:-0.085}" \
    FIG6_SMAX="${FIG6_SMAX:-0.120}" \
    FIG6_DS="${FIG6_DS:-0.001}" \
    FIG6_RELAX="${FIG6_RELAX:-80.0}" \
    FIG6_FREE="${FIG6_FREE:-80.0}" \
    FIG6_EPS="${FIG6_EPS:-0.001}"
