# NeuroThermo publication release v1.2

This directory is the self-contained data, code, and result release for the
WT--SCA3 Purkinje-cell state-space manuscript. It contains the processed data
used in the article, executable analysis code, publication figures, and frozen
validation records.

## Reproduce the published figures

From this directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
./code/reproduce_figures.sh
```

The generated files are written to `results/figures/`. Figures 1--3 and S1 can
also be regenerated with R:

```bash
Rscript code/figures/R/make_all_figures.R
```

## Validate the release

```bash
python code/validate_release.py
./code/run_smoke_tests.sh
```

The first command checks file integrity and all central numerical claims. The
second runs unit and reduced-data smoke tests for both simulation pipelines.

## Re-run the full simulations

```bash
./code/run_full_analyses.sh kl
./code/run_full_analyses.sh nonequilibrium
```

Both commands resolve inputs from `data/inputs/`; no external or
machine-specific data path is used. These are computationally intensive runs.

## Directory map

- `data/inputs/`: analysis-ready inputs for both full simulation pipelines.
- `data/figure_source/`: numerical source data for publication figures and
  sensitivity tables.
- `data/kl_convergence_v1_0_1/`: full-coverage KL outputs.
- `data/nonequilibrium_geometry_v1_0_1/`: article-facing nonequilibrium outputs.
- `code/pipelines/`: complete versioned simulation packages and tests.
- `code/figures/`: Python and R figure-generation code.
- `results/figures/`: publication-ready PDF and PNG figures.
- `results/validation/`: frozen verdicts, numerical QC, and run fingerprints.
- `docs/`: reproducibility notes and result-to-code traceability.

The raw historical ABF archive is not duplicated here. The release contains
the de-identified, analysis-ready tables required for the reported analyses.
See `data/README.md` for cohort and provenance details.
