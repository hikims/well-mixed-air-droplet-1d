# GitHub upload checklist

## Before committing

- [ ] Remove compiled binaries from Git tracking.
- [ ] Remove full simulation folders unless intentionally archived.
- [ ] Keep only source code, plotting scripts, small processed CSV files, and final figure outputs.
- [ ] Check that `README.md` reproduces the actual commands.
- [ ] Check that `Data availability` in the manuscript points to this repository.
- [ ] Ask coauthor to confirm license.
- [ ] Run `make figures` from a fresh clone or copied folder.

## Suggested repository name

```text
well-mixed-air-droplet-1d
```

## Manuscript data availability text

```latex
The simulation code, plotting scripts, and processed data needed to reproduce the figures in this paper are openly available at \url{https://github.com/hikims/well-mixed-air-droplet-1d}.
```

## First GitHub push

```bash
git init
git add README.md LICENSE CITATION.cff .gitignore environment.yml requirements.txt Makefile
git add src include scripts docs
git add data/processed figures/main figures/supplement
git commit -m "Initial reproducibility release for well-mixed-air droplet model"

git branch -M main
git remote add origin git@github.com:hikims/well-mixed-air-droplet-1d.git
git push -u origin main
```

If using HTTPS:

```bash
git remote add origin https://github.com/hikims/well-mixed-air-droplet-1d.git
```
