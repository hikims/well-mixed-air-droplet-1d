#!/usr/bin/env bash
set -euo pipefail

# Grid-aligned nonlinear validation near the first transition.
# Reproduces manuscript Fig. 8 and Supplementary Figs. S4-S6.
#
# Protocol (manuscript Section 4.7):
#   Nx = 401, dt = 6.25e-5, T_relax = 100, T_free = 400,
#   position perturbations eps = 1e-3 and 1e-5,
#   late-time diagnostics over the final 25% of the free interval.
#
# Requires the Hopf reference outputs first:
#   make hopf
# (results/critical_points.csv and results/time_step_stability.csv)
#
# Usage from the repository root:
#   bash scripts/run_grid_aligned_validation.sh
# Override any parameter via environment variables, e.g.
#   THREADS=12 bash scripts/run_grid_aligned_validation.sh

OUT="${OUT:-data/generated}"
BIN="${BIN:-bin/grid_aligned_validation}"
SMIN="${SMIN:-0.081750}"
SMAX="${SMAX:-0.082750}"
DS0="${DS0:-0.000125}"
TRELAX="${TRELAX:-100.0}"
TFREE="${TFREE:-400.0}"
DT="${DT:-6.25e-5}"
SAVE_EVERY="${SAVE_EVERY:-20}"
THREADS="${THREADS:-8}"

HOPF_CSV="${HOPF_CSV:-results/critical_points.csv}"
IMEX_CSV="${IMEX_CSV:-results/time_step_stability.csv}"
FIG_DIR="${FIG_DIR:-figures/supplement/grid_aligned}"

if [[ ! -f "${HOPF_CSV}" || ! -f "${IMEX_CSV}" ]]; then
    echo "ERROR: Hopf reference files not found:" >&2
    echo "  ${HOPF_CSV}" >&2
    echo "  ${IMEX_CSV}" >&2
    echo "Run 'make hopf' first." >&2
    exit 1
fi

make "${BIN}"
mkdir -p "${OUT}"

declare -A EPS_VALUES=( ["eps1e-3"]="1.0e-3" ["eps1e-5"]="1.0e-5" )

for label in "eps1e-3" "eps1e-5"; do
    eps="${EPS_VALUES[$label]}"
    out_dir="${OUT}/grid_aligned_${label}"
    echo "=== grid-aligned sweep ${label} (eps=${eps}) -> ${out_dir} ==="
    "${BIN}" "${out_dir}" "${SMIN}" "${SMAX}" "${DS0}" \
        "${TRELAX}" "${TFREE}" "${DT}" "${SAVE_EVERY}" "${THREADS}" "${eps}"
done

python3 scripts/analyze_grid_aligned_validation.py \
    --scan "eps1e-3=${OUT}/grid_aligned_eps1e-3" \
    --scan "eps1e-5=${OUT}/grid_aligned_eps1e-5" \
    --output-dir "${FIG_DIR}" \
    --hopf-csv "${HOPF_CSV}" \
    --imex-csv "${IMEX_CSV}"

echo "Done. Figures written to ${FIG_DIR}"
