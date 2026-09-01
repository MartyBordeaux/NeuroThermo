# NeuroThermo publication release v1.2

This directory is the reproducibility release for the manuscript *Multiscale dynamical staging in WT--SCA3 Purkinje-cell models*. It contains the 50 raw current-clamp ABF recordings, QC/provenance material, frozen manual selections, exact analysis source, computationally heavy frozen checkpoints/results, publication-facing source tables, and canonical figure renderers.

All publication configs use repository-relative paths. The tested Python environment is Python 3.9.25.

## 1. Create the publication environment

From `publication/neurothermo_plos_v1_2/`:

```bash
python3.9 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r environment/requirements-python39-historical.txt
```

The portable publication lock pins NumPy 1.25.2, pandas 2.2.3, SciPy 1.13.1, scikit-learn 1.6.1, pyABF 2.3.8, numba 0.60.0 and llvmlite 0.43.0. NumPy 1.25.2 is intentional: cross-runner diagnostics found CPU-dependent candidate identities under NumPy 1.26.4, whereas 1.25.2 reproduced all 6800 frozen identities across the tested hosted runners and SciPy versions.

The R plotting scripts are retained as optional alternative renderers. R is not required for the canonical clean-clone reproduction path; publication figures are reproduced with Python/matplotlib under the pinned environment above.

## 2. Run the release integrity gate

```bash
bash code/run_full_analyses.sh preflight
bash code/run_full_analyses.sh raw-integrity
bash code/run_full_analyses.sh qc-tests
```

These commands verify the frozen calibration archive, full imported v3.9/dynamic inputs, 50 raw ABFs and their SHA-256 manifest, exact transition frozen inputs, and cross-stage transition-result integrity.

## 3. Recompute raw ABF -> fixed QC

```bash
bash code/run_full_analyses.sh qc-recompute
```

Under the pinned environment this reproduces the publication QC layer: 6800 candidate events and 6039 final `fixed_qc_detected` spikes. The final fixed-QC implementation includes the frozen manual overrides rather than substituting a new spike detector or threshold.

## 4. Restricted four-parameter HR fit

First validate the exact frozen fit inputs:

```bash
bash code/run_full_analyses.sh cellfit-validate
```

To perform the full expensive v3.9 refit and subsequent practical-identifiability calculation:

```bash
bash code/run_full_analyses.sh cellfit
bash code/run_full_analyses.sh cellfit-identify
```

The model fits `b`, `r`, `s`, and `kappa_I`. Spike times are fitted with the Victor--Purpura objective; non-spiking sweeps enter only through the binary rheobase bracket. Per-sweep latency alignment is an exact additive first-spike shift. The last spike is not anchored and time is not rescaled between first and last spikes.

The full fit is intentionally not run in routine hosted CI because it is computationally expensive. Its exact frozen result tree is included under `data/v3_9_results_full/` so downstream analyses can be reproduced without rerunning the optimizer.

## 5. Post-fit and endpoint analyses

The publication replay path uses the exact frozen v3.9 results and the published frozen inputs for each downstream stage:

```bash
bash code/run_full_analyses.sh characterization
bash code/run_full_analyses.sh dynamic-validate
bash code/run_full_analyses.sh dynamic
bash code/run_full_analyses.sh endpoint-validate
bash code/run_full_analyses.sh endpoint
```

`characterization` deterministically creates the historical XLSX compatibility view from `data/animal_id_recovery/accepted_cohort.csv`; the CSV remains the scientific source of truth.

## 6. Transition analysis

Lightweight validation/reprojection sequence:

```bash
bash code/run_full_analyses.sh transition-frozen
bash code/run_full_analyses.sh transition-v1-0-validate
bash code/run_full_analyses.sh transition-v1-1-validate
bash code/run_full_analyses.sh transition-v1-1
bash code/run_full_analyses.sh prepare-transition-v1-2
bash code/run_full_analyses.sh transition-v1-2-validate
bash code/run_full_analyses.sh transition-v1-2-1-validate-frozen
bash code/run_full_analyses.sh transition-v1-3-validate-frozen
bash code/run_full_analyses.sh transition-integrity
```

This sequence is clean-clone tested. The v1.1 reprojection performs zero new HR simulations, and the generated v1.1 -> v1.2 staging inputs are checked against their historical references.

For a full expensive transition recomputation, use:

```bash
bash code/run_full_analyses.sh transition-v1-0
bash code/run_full_analyses.sh transition-v1-1
bash code/run_full_analyses.sh prepare-transition-v1-2
bash code/run_full_analyses.sh transition-v1-2
bash code/run_full_analyses.sh transition-v1-2-1
bash code/run_full_analyses.sh transition-v1-3
```

The v1.2 and v1.3 stages are checkpoint/resume calculations. The release also includes the exact historical checkpoint/result trees, allowing validation and downstream reproduction without recomputing approximately 0.95 million v1.2 and 1.90 million new v1.3 HR states in routine CI.

## 7. KL and nonequilibrium/Fisher analyses

```bash
bash code/run_full_analyses.sh kl
bash code/run_full_analyses.sh nonequilibrium
```

The KL package contains the full-state and marginal convergence analysis used for the KL-before-firing result. The nonequilibrium-geometry package contains the Fisher/time-reversal/stationarity analysis used by the manuscript; a second publication-critical Fisher package is not required.

## 8. Rebuild publication source data and figures

The single figure-reproduction command is:

```bash
bash code/reproduce_figures.sh
```

Before rendering, this command runs `code/assemble_figure_source.py`, which reconstructs all eight Fig.1--3 source tables actually consumed by the canonical renderer and checks them against the immutable references in `data/figure_source/` with a numerical tolerance of 1e-12. Recomputed copies are written to `results/recomputed/figure_source/`; the frozen references are never overwritten.

The canonical renderer writes Figures 1--5 and S1--S3 to `results/figures/`, each in PDF and PNG form. The verifier requires exactly these 16 freshly generated files and rejects missing or extra renders.

Publication-facing capacitance normalization is explicit in the source-data assembler: capacitance and `kappa_I` are doubled, current density `J=I/Cm` is halved, and current in pA is unchanged. For the four WT recordings whose animal identity is not recoverable, Fig.1 recording-date labels are derived directly from the raw ABF headers; they are not treated as recovered animal IDs.

## Two reproducibility levels

The release deliberately separates two tasks:

1. **Clean-clone publication replay.** Exact frozen checkpoints/results are supplied for computationally expensive stages. All source-data assembly, integrity gates, lightweight reprojections, and figure rendering can therefore be checked without days of optimization/simulation.
2. **Full expensive recomputation.** Explicit `cellfit`, `cellfit-identify`, `transition-v1-0`, `transition-v1-2`, and `transition-v1-3` commands are provided with the exact configs, seeds, bounds and checkpoint/resume behavior used by the publication analysis.

A frozen intermediate is not silently substituted: its origin and SHA-256 are recorded in the release provenance files and workflow inventory.

## Directory map

- `data/raw/WT`, `data/raw/SCA3`: 50 raw ABF recordings.
- `data/qc/`: QC/provenance workbook and its fingerprint.
- `data/calibration_bundle/`: exact server calibration archive.
- `data/calibration/`: deterministic extracted frozen fit/QC inputs plus `CALIBRATION_PROVENANCE.tsv`.
- `data/v3_9_results_full/`: full exact v3.9 frozen result tree for downstream replay.
- `data/transition_v1_*_frozen/`: exact transition frozen inputs.
- `data/transition_v1_*_results/`: exact transition result/checkpoint trees.
- `data/kl_convergence_v1_0_1/`: publication KL outputs.
- `data/nonequilibrium_geometry_v1_0_1/`: publication nonequilibrium/Fisher outputs.
- `data/figure_source/`: immutable publication-facing source tables.
- `code/pipelines/`: final versioned analysis packages.
- `code/assemble_figure_source.py`: deterministic Fig.1--3 source-data assembly and frozen-equivalence gate.
- `code/run_full_analyses.sh`: publication master command dispatcher.
- `code/figures/python/`: canonical figure-generation code.
- `code/figures/R/`: optional alternative R renderers.
- `results/recomputed/`: outputs generated by clean-clone replay/recomputation.
- `results/figures/`: publication-ready figures.
- `docs/PUBLICATION_WORKFLOW_INVENTORY.tsv`: stage-by-stage computational inventory.
- `docs/RESULT_TO_CODE_MAP.md`: article result -> upstream pipeline -> source data -> renderer mapping.

The constructed path coordinate `p` is a coordinate through model-state space. It is not interpreted as disease time or as evidence of an irreversible biological trajectory.
