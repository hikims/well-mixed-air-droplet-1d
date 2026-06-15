# Reproducibility guide

## 1. Build

```bash
make all
```

## 2. Main manuscript figures

```bash
make figures
```

Expected outputs:

```text
figures/main/figure1_phase_diagram.png
figures/main/figure1_phase_diagram.eps
figures/main/figure2_representative_trajectories_compact.png
figures/main/figure2_representative_trajectories_compact.eps
...
```

## 3. Supplementary checks

```bash
make supplement
```

Expected checks:

- grid refinement
- time-step refinement
- long-time integration near the first transition
- perturbation-amplitude check
- cover-leakage sensitivity

## 4. Notes on large data

Do not commit full raw output directories unless the repository is intended as an archival data repository.  For GitHub, keep only lightweight processed CSV files needed for plotting, or use GitHub Releases / Zenodo for larger data.
