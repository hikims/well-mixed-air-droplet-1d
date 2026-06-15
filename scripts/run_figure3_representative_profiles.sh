#!/usr/bin/env bash
set -e

python3 scripts/plot_figure3_representative_profiles.py "${1:-figure2_data}" "${2:-figure}" "${3:-1.0}"
