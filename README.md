# Well-mixed-air reduction for self-propelled droplets

This repository contains simulation and plotting code for the manuscript

**A well-mixed-air reduction to a one-dimensional reaction--diffusion particle model for self-propelled oil droplets**

by Minsoo Kim and Masaharu Nagayama.

The code reproduces the single-droplet phase diagram, representative trajectories, chemical profiles, surface-tension imbalance diagnostics, cover-protocol simulation, the direct odd-sector Hopf-root calculation H(lambda)=1, the IMEX time-step stability analysis, the grid-aligned nonlinear validation near onset, and the supplementary robustness checks.

## Repository status

This repository is intended for reproducibility of the manuscript figures. Large generated data files are not tracked by default. Regenerate them using the Makefile targets below, or place archived data under `data/processed/` if you want a lightweight reproducibility snapshot.

## Quick start

```bash
git clone https://github.com/hikims/well-mixed-air-droplet-1d.git
cd well-mixed-air-droplet-1d

# Python environment
conda env create -f environment.yml
conda activate droplet1d

# Build C++ binaries
make all

# Check scripts and file-name consistency
make check-python
make check-names

# Run simulations and generate all main manuscript figures
make figures

# Run supplementary robustness checks
make supplement
```

## Main Makefile targets

```bash
make show          # show detected parameters and build settings
make all           # build all main C++ binaries
make figures       # generate all main manuscript figures
make supplement    # run supplementary robustness checks
make check-python  # check Python script syntax
make check-names   # check for obsolete file-name references
make clean         # remove compiled binaries
make clean-data    # remove generated data and logs
make distclean     # remove binaries, generated data, logs, and generated figures
```

Individual main figures can be regenerated with:

```bash
make figure1       # phase diagram                         (manuscript Fig. 2)
make figure2       # representative trajectories            (manuscript Fig. 3)
make figure3       # representative chemical profiles       (manuscript Fig. 4)
make figures4_5    # surface-tension imbalance              (manuscript Figs. 5-6)
make figure6       # cover protocol                         (manuscript Fig. 7, Fig. S1)
make hopf          # odd-sector Hopf roots + IMEX stability (manuscript Fig. 9, S2-S3, Tables S1-S2)
make grid_aligned  # grid-aligned validation, after 'make hopf' (manuscript Fig. 8, S4-S6)
```

Internal `figureN` target numbers predate the final manuscript numbering
(the manuscript's Fig. 1 is a schematic with no associated code); the
manuscript mapping in parentheses above is authoritative.

Long runs can also be launched in the background:

```bash
make figure1_bg
```

## Main outputs

Generated manuscript figures are written to

```text
figures/main/
figures/supplement/
```

Simulation outputs are written to

```text
data/generated/
data/processed/
```

The final manuscript uses the following figure files:

| Manuscript item | File |
| --- | --- |
| Fig. 2 (phase diagram) | `figures/main/figure1_phase_diagram.eps` |
| Fig. 3 (trajectories) | `figures/main/figure2_representative_trajectories_compact.eps` |
| Fig. 4 (chemical profiles) | `figures/main/figure3_representative_profiles_compact.eps`, `figures/main/figure3_sds_profiles_zoom.eps` |
| Fig. 5 (oscillatory imbalance) | `figures/main/figure4_oscillatory_surface_tension_detail.eps` |
| Fig. 6 (rotating imbalance) | `figures/main/figure4_rotating_surface_tension_detail.eps` |
| Fig. 7 (cover protocol) | `figures/main/figure5_cover_protocol_clean.eps` |
| Fig. 8, Figs. S4-S6 (grid-aligned onset) | `figures/supplement/grid_aligned/` |
| Fig. 9, Figs. S2-S3 (Hopf roots, time-step shift) | `results/figures/` |
| Fig. S1 (cover-stage profiles) | cover-protocol stage profiles from `make figure6` |
| Figs. S7-S8 (grid refinement, cover leakage) | `figures/supplement/` via `make supplement` |
| Tables S1-S2 | `results/critical_points.csv`, `results/time_step_stability.csv` |

## Repository layout

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── requirements.txt
├── Makefile
├── src/
│   ├── model_lu_tbb.hpp
│   ├── params1.h
│   ├── figure1_phase_diagram.cpp
│   ├── figure2_trajectories.cpp
│   ├── figure4_surface_tension.cpp
│   ├── figure5_cover_protocol.cpp
│   ├── grid_aligned_validation.cpp
│   ├── hopf_solver.py
│   ├── imex_stability.py
│   ├── supplement_defense_batch.cpp
│   ├── supplement_oscillatory_cycle.cpp
│   └── legacy/
│       ├── figure6_perturbation_activity.cpp   # superseded half-grid perturbation test
│       ├── amplitude_lu.cpp
│       ├── close_lu.cpp
│       ├── open_lu.cpp
│       ├── reopen_lu.cpp
│       └── rest_perturb_lu.cpp
├── scripts/
│   ├── patch_params.py
│   ├── plot_figure1_phase_diagram.py
│   ├── plot_figure2_representative_trajectories.py
│   ├── plot_figure3_representative_profiles.py
│   ├── plot_figure4_5_surface_tension.py
│   ├── plot_figure6_cover_protocol.py
│   ├── analyze_grid_aligned_validation.py
│   ├── plot_supplement_robustness_checks.py
│   ├── run_figure1_phase_diagram.sh
│   ├── run_figure2_representative_trajectories.sh
│   ├── run_figure3_representative_profiles.sh
│   ├── run_figure4_5_surface_tension.sh
│   ├── run_figure6_cover_protocol.sh
│   ├── run_grid_aligned_validation.sh
│   ├── run_hopf_scan.py
│   ├── run_time_step_stability.py
│   ├── plot_hopf_results.py
│   └── run_supplement_robustness_checks.sh
├── data/
│   ├── README.md
│   ├── generated/
│   └── processed/
├── figures/
│   ├── README.md
│   ├── main/
│   └── supplement/
└── docs/
    ├── DATA_DICTIONARY.md
    ├── GITHUB_UPLOAD_CHECKLIST.md
    ├── REPO_STRUCTURE.md
    └── REPRODUCIBILITY.md
```

## Reproducibility notes

The default numerical parameters follow the manuscript parameter table:

* `Nx = 401`
* `dt = 5.0e-4`
* `L = 16`
* `r = 1`
* open-cover leakage `kw_open = 10000`
* closed-cover leakage `kw_closed = 0`

The main parameter header is

```text
src/params1.h
```

The supplementary robustness checks include grid refinement, time-step refinement, long-time checks near the first transition, perturbation-amplitude checks, and cover-leakage sensitivity.

## Notes on legacy files

The `src/legacy/` directory contains earlier single-purpose drivers used during model development. These files are not required to reproduce the final manuscript figures, but are retained for transparency.

## License

Released under the MIT License; see `LICENSE`.
