#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# One-command defense simulations for the Physica D manuscript
# ============================================================
# Run from the project root:
#   bash scripts/run_defense_batch.sh
#
# Useful environment variables:
#   OUT=defense_data THREADS=8 USE_TBB=1 bash scripts/run_defense_batch.sh
#   RUN_LONG=0 bash scripts/run_defense_batch.sh       # skip long-time cases
#   RUN_TRANSITION=0 bash scripts/run_defense_batch.sh # skip fine transition scan
#
# This script edits the parameter header to change Nx, rebuilds the defense binary,
# runs the checks, and restores the original parameter header when finished.
# It auto-detects src/params.h or src/params1.h. You can override with:
#   PARAMS_FILE=src/params1.h bash scripts/run_defense_batch.sh

OUT="${OUT:-defense_data}"
THREADS="${THREADS:-8}"
USE_TBB="${USE_TBB:-1}"
KW_OPEN_DEFAULT="${KW_OPEN_DEFAULT:-10000}"
KW_CLOSE_DEFAULT="${KW_CLOSE_DEFAULT:-0}"

RUN_GRID="${RUN_GRID:-1}"
RUN_TIME="${RUN_TIME:-1}"
RUN_TRANSITION="${RUN_TRANSITION:-1}"
RUN_LONG="${RUN_LONG:-1}"
RUN_PERTURB="${RUN_PERTURB:-1}"
RUN_COVER="${RUN_COVER:-1}"

mkdir -p "${OUT}" logs

# ------------------------------------------------------------
# Parameter header auto-detection
# ------------------------------------------------------------
PARAMS_FILE="${PARAMS_FILE:-}"
if [[ -z "${PARAMS_FILE}" ]]; then
    if [[ -f src/params.h ]]; then
        PARAMS_FILE="src/params.h"
    elif [[ -f src/params1.h ]]; then
        PARAMS_FILE="src/params1.h"
    else
        echo "ERROR: run this script from the project root; neither src/params.h nor src/params1.h was found." >&2
        echo "       Existing parameter headers:" >&2
        ls src/params*.h 2>/dev/null >&2 || true
        exit 1
    fi
fi

if [[ ! -f "${PARAMS_FILE}" ]]; then
    echo "ERROR: PARAMS_FILE=${PARAMS_FILE} does not exist." >&2
    exit 1
fi

echo "Using parameter header: ${PARAMS_FILE}"

PARAM_BAK="${OUT}/$(basename "${PARAMS_FILE}").backup.$(date +%Y%m%d_%H%M%S)"
cp "${PARAMS_FILE}" "${PARAM_BAK}"
restore_params() {
    cp "${PARAM_BAK}" "${PARAMS_FILE}"
}
trap restore_params EXIT

build_for_nx() {
    local nx="$1"
    echo "=== Building defense binary for Nx=${nx} ==="
    python3 scripts/patch_params.py --params "${PARAMS_FILE}" --nx "${nx}"
    make clean >/dev/null 2>&1 || true
    make defense_batch USE_TBB="${USE_TBB}" THREADS="${THREADS}"
}

run_open() {
    local suite="$1"; shift
    local dt="$1"; shift
    local tend="$1"; shift
    local late="$1"; shift
    local kw="$1"; shift
    echo "--- open suite=${suite} dt=${dt} tend=${tend} ---"
    TBB_NUM_THREADS="${THREADS}" ./bin/supplement_defense_batch open "${OUT}" "${suite}" "${dt}" "${tend}" "${late}" "${kw}" "$@"
}

run_perturb() {
    local suite="$1"; shift
    local dt="$1"; shift
    local trelax="$1"; shift
    local tfree="$1"; shift
    local late="$1"; shift
    local kw="$1"; shift
    local eps="$1"; shift
    echo "--- perturb suite=${suite} eps=${eps} ---"
    TBB_NUM_THREADS="${THREADS}" ./bin/supplement_defense_batch perturb "${OUT}" "${suite}" "${dt}" "${trelax}" "${tfree}" "${late}" "${kw}" "${eps}" "$@"
}

run_cover() {
    local suite="$1"; shift
    local dt="$1"; shift
    local topen="$1"; shift
    local tclosed="$1"; shift
    local treopen="$1"; shift
    local kwopen="$1"; shift
    local kwclose="$1"; shift
    local s0="$1"; shift
    local seed="$1"; shift
    echo "--- cover suite=${suite} kw_open=${kwopen} ---"
    TBB_NUM_THREADS="${THREADS}" ./bin/supplement_defense_batch cover "${OUT}" "${suite}" "${dt}" "${topen}" "${tclosed}" "${treopen}" "${kwopen}" "${kwclose}" "${s0}" "${seed}"
}

# ------------------------------------------------------------
# 1. Grid refinement: representative states and near-boundary cases
# ------------------------------------------------------------
if [[ "${RUN_GRID}" == "1" ]]; then
    for nx in 401 801; do
        build_for_nx "${nx}"
        run_open "grid_Nx${nx}_representative" 0.0005 220 40 "${KW_OPEN_DEFAULT}" 0.080 0.125 0.180
        run_open "grid_Nx${nx}_boundary1"      0.0005 260 60 "${KW_OPEN_DEFAULT}" 0.098 0.100 0.102
        run_open "grid_Nx${nx}_boundary2"      0.0005 260 60 "${KW_OPEN_DEFAULT}" 0.148 0.152 0.156
    done
fi

# ------------------------------------------------------------
# 2. Time-step refinement at baseline grid
# ------------------------------------------------------------
if [[ "${RUN_TIME}" == "1" ]]; then
    build_for_nx 401

    # Use representative states away from transition boundaries.
    # dt = 1e-3 was too coarse, so we test finer time steps.
    run_open "dt_5e-4"     0.0005   220 40 "${KW_OPEN_DEFAULT}" 0.080 0.125 0.180
    run_open "dt_2p5e-4"   0.00025  220 40 "${KW_OPEN_DEFAULT}" 0.080 0.125 0.180
    run_open "dt_1p25e-4"  0.000125 220 40 "${KW_OPEN_DEFAULT}" 0.080 0.125 0.180
fi

# ------------------------------------------------------------
# 3. Fine transition scans
# ------------------------------------------------------------
if [[ "${RUN_TRANSITION}" == "1" ]]; then
    build_for_nx 401
    run_open "transition_first"  0.0005 320 80 "${KW_OPEN_DEFAULT}" \
        0.095 0.096 0.097 0.098 0.099 0.100 0.101 0.102 0.103 0.104 0.105
    run_open "transition_second" 0.0005 320 80 "${KW_OPEN_DEFAULT}" \
        0.145 0.147 0.149 0.151 0.153 0.155 0.157 0.159
fi

# ------------------------------------------------------------
# 4. Long-time checks near the first transition
# ------------------------------------------------------------
if [[ "${RUN_LONG}" == "1" ]]; then
    build_for_nx 401
    run_open "longtime_500"  0.0005 500  100 "${KW_OPEN_DEFAULT}" 0.098 0.099 0.100 0.101 0.102
    run_open "longtime_1000" 0.0005 1000 200 "${KW_OPEN_DEFAULT}" 0.098 0.099 0.100 0.101 0.102
fi

# ------------------------------------------------------------
# 5. Perturbation amplitude and sign checks
# ------------------------------------------------------------
if [[ "${RUN_PERTURB}" == "1" ]]; then
    build_for_nx 401
    for eps in 0.005 0.010 0.020 -0.010; do
        label="eps_${eps}"
        label="${label//- /m}"
        run_perturb "perturb_${label}" 0.0005 80 80 20 "${KW_OPEN_DEFAULT}" "${eps}" 0.095 0.100 0.105
    done
fi

# ------------------------------------------------------------
# 6. Cover leakage sensitivity
# ------------------------------------------------------------
if [[ "${RUN_COVER}" == "1" ]]; then
    build_for_nx 401
    run_cover "cover_kw_1e3"  0.0005 40 40 40 1000   "${KW_CLOSE_DEFAULT}" 0.125 3101
    run_cover "cover_kw_1e4"  0.0005 40 40 40 10000  "${KW_CLOSE_DEFAULT}" 0.125 3102
    run_cover "cover_kw_1e5"  0.0005 40 40 40 100000 "${KW_CLOSE_DEFAULT}" 0.125 3103
fi

# Restore original params before plotting.
restore_params
trap - EXIT

python3 scripts/plot_supplement_robustness_checks.py "${OUT}"

echo "=== Defense batch complete ==="
echo "Results: ${OUT}"
echo "Summary figures: ${OUT}/figures"
echo "LaTeX snippets: ${OUT}/latex_snippets.tex"
