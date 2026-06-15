# Proposed repository structure

Use this structure as the final public layout.

```text
well-mixed-air-droplet-1d/
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── requirements.txt
├── Makefile
├── include/
│   └── params.h
├── src/
│   ├── common/
│   ├── figure1_phase_diagram.cpp
│   ├── figure2_trajectories.cpp
│   ├── figure3_profiles.cpp
│   ├── figure4_surface_tension.cpp
│   ├── figure5_cover_protocol.cpp
│   ├── figure6_perturbation_activity.cpp
│   └── defense_batch.cpp
├── scripts/
│   ├── run_all_main.sh
│   ├── run_all_supplement.sh
│   ├── plot_figure1_phase_diagram.py
│   ├── plot_figure2_representative_trajectories.py
│   ├── plot_figure3_representative_profiles.py
│   ├── plot_figure4_5_surface_tension.py
│   ├── plot_figure6_cover_protocol.py
│   ├── plot_figure6_perturbation_activity.py
│   └── plot_supplement_robustness_checks.py
├── data/
│   ├── processed/
│   └── README.md
├── figures/
│   ├── main/
│   └── supplement/
└── docs/
    ├── REPRODUCIBILITY.md
    ├── DATA_DICTIONARY.md
    └── GITHUB_UPLOAD_CHECKLIST.md
```
