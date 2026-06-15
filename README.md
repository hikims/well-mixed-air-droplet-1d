````markdown
# Well-mixed-air reduction for self-propelled droplets

This repository contains simulation and plotting code for the manuscript

**A well-mixed-air reduction of a one-dimensional reaction--diffusion particle model for self-propelled oil droplets**

by Minsoo Kim and Masaharu Nagayama.

The code reproduces the single-droplet phase diagram, representative trajectories, chemical profiles, surface-tension imbalance diagnostics, cover-protocol simulation, perturbation-activity test, and supplementary robustness checks.

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
````

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
make figure1       # phase diagram
make figure2       # representative trajectories
make figure3       # representative chemical profiles
make figures4_5    # oscillatory and rotating surface-tension imbalance
make figure6       # cover protocol
make figure7       # perturbation activity near transition
```

Long runs can also be launched in the background:

```bash
make figure1_bg
make figure7_bg
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

The manuscript currently uses the following main figure files:

| Figure | File                                                                                                     |
| ------ | -------------------------------------------------------------------------------------------------------- |
| Fig. 1 | `figures/main/figure1_phase_diagram.eps`                                                                 |
| Fig. 2 | `figures/main/figure2_representative_trajectories_compact.eps`                                           |
| Fig. 3 | `figures/main/figure3_representative_profiles_compact.eps`, `figures/main/figure3_sds_profiles_zoom.eps` |
| Fig. 4 | `figures/main/figure4_oscillatory_surface_tension_detail.eps`                                            |
| Fig. 5 | `figures/main/figure4_rotating_surface_tension_detail.eps`                                               |
| Fig. 6 | `figures/main/figure5_cover_protocol_clean.eps`                                                          |
| Fig. 7 | `figures/main/figure6_activity_onset.eps`                                                                |

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
│   ├── figure6_perturbation_activity.cpp
│   ├── supplement_defense_batch.cpp
│   ├── supplement_oscillatory_cycle.cpp
│   └── legacy/
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
│   ├── plot_figure7_perturbation_activity.py
│   ├── plot_supplement_robustness_checks.py
│   ├── run_figure1_phase_diagram.sh
│   ├── run_figure2_representative_trajectories.sh
│   ├── run_figure3_representative_profiles.sh
│   ├── run_figure4_5_surface_tension.sh
│   ├── run_figure6_cover_protocol.sh
│   ├── run_figure7_perturbation_activity.sh
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

## Data availability statement for manuscript

After creating the GitHub repository, replace the manuscript data availability statement with:

```latex
The simulation code, plotting scripts, and processed data needed to reproduce the figures in this paper are openly available at \url{https://github.com/hikims/well-mixed-air-droplet-1d}.
```

Replace the URL if you choose a different repository name.

## License

See `LICENSE`. Confirm the license with all coauthors before public release.

```
```
