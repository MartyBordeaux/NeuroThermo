# NeuroThermo

Frozen analysis repository for Hindmarsh–Rose ensemble fitting, uncertainty-aware dynamical characterization, and WT→SCA3 transition staging of Purkinje-cell current-clamp recordings.

## Current frozen stage

The active analysis is now:

1. **cell-fit v3.9** — final four-parameter HR cell fits;
2. **post-fit characterization v1.0** — animal-aware parameter and phenotype characterization;
3. **dynamic characterization v2.1** — experiment-support-restricted suprathreshold dynamics and near-optimal-solution robustness;
4. **endpoint ensemble v1.0.1** — uncertainty-aware WT and SCA3 endpoint clouds;
5. **transition ensemble v1.1** — corrected latency-invariant staging geometry;
6. **transition ensemble v1.2.1** — scenario-first intrinsic × drive map;
7. **transition ensemble v1.3.1** — factorized drive decomposition into model input scaling `kappa_I`, experimental current protocol `J`, and their interaction. v1.3.1 is a figure-only correction of v1.3.0; numerical v1.3 results are unchanged.

### Frozen cell-fit model

- Fitted cell-level HR coordinates: `b`, `r`, `s`, `kappa_I`.
- Fixed HR constants: `a=1`, `c=1`, `d=5`, `x_R=-1.6`.
- Primary cohort: **18 multi-sweep cells (12 WT, 6 SCA3)**.
- Secondary descriptive set: 2 accepted single-sweep cells.
- Exact additive first-spike alignment per spiking sweep; no time rescaling and no last-spike anchoring.
- Binary rheobase bracket retained as the absolute excitability constraint.
- Final common bounds: `b=[0.5,7]`, `r=[1e-4,0.1]`, `s=[0.05,15]`, `kappa_I=[2e-4,2]`.

### Frozen dynamical phenotype

Dynamic v2.1 restricts primary comparisons to experimentally supported suprathreshold regions. The shared progression coordinate is

`q = (J - J_rheo) / (J_max,obs - J_rheo)`.

No extrapolation beyond the observed spiking-current range is allowed in the primary analysis. The strongest common support is at `q=0.75` (18/18 cells), with `q=0.50` available in 17/18 cells.

The primary robust phenotype is higher capacitance-normalized rheobase in SCA3 together with slower experiment-supported suprathreshold dynamics.

### Frozen transition geometry

The primary transition projection is based on `(log10 J_rheo, log10 mean_ISI_q75)`. In the core-secure endpoint set:

- WT-exit boundary: `A_ISI = 0.135829`;
- balance: `A_ISI = 0.5`;
- SCA3-entry boundary: `A_ISI = 0.797856`.

For the 32 core-secure WT×SCA3 biological pairs, coupled-path median staging is approximately `p=0.398` (WT-exit), `0.676` (balance), and `0.837` (SCA3-entry). Early versus late drive timing shifts the transition substantially.

Scenario-first 2D analysis shows stage-dependent control: intrinsic and drive contributions are comparable near WT-exit, nearly balanced around the transition midpoint, and the drive contribution increases toward SCA3-entry.

### Drive decomposition

The combined drive used in v1.2 is decomposed in v1.3 into:

- `kappa_I`: fitted HR input-scaling coordinate;
- `J`: experimenter-imposed current protocol;
- a non-additive interaction term.

The SCA3-directed effect is primarily associated with `kappa_I`; changing `J` alone does not reproduce the same transition and generally opposes part of the `kappa_I` effect. `J` is therefore treated as protocol sensitivity, not a disease parameter. Raw `kappa_I` is also not interpreted as a direct biophysical conductance because of its known association with capacitance.

## Repository layout

```text
docs/                                  frozen methods, scientific status, provenance
results/v3_9/                           final cell parameters
results/characterization_v1/            animal-aware characterization
results/dynamic_v2_1/                   support-restricted dynamic characterization
results/endpoint_ensemble_v1_0_1/       transition-ready endpoint ensemble
results/transition_v1_1/                primary corrected staging
results/transition_v1_2_1/              scenario-first 2D intrinsic × drive analysis
results/transition_v1_3/                drive decomposition and corrected figures
```

Earlier working pipelines remain recoverable from Git history. The current branch is intentionally a compact frozen-analysis snapshot rather than a chronological dump of intermediate calculations.

See `docs/SCIENTIFIC_STATUS.md`, `docs/TRANSITION_V1_1_TO_V1_3.md`, `docs/DYNAMIC_V2_1.md`, `docs/ENDPOINT_ENSEMBLE_V1_0_1.md`, `docs/METHODS_V3_9.md`, and `docs/ANIMAL_ID_PROVENANCE.md`.

## Parallel raw-data phenotype track

A complementary model-free analysis derives cell-level phenotypes directly from the same current-clamp recordings. Its frozen stage is **raw-data phenotype v0.7.3**, which identifies a short-timescale 4–8 ms ordinal temporal-order deficit in SCA3 that survives an exact-spectrum rank-Gaussian Fourier null.

This track retains all 13 WT and 7 SCA3 cells and uses the recorded cell as the independent unit. Animal-level inference is absent. Current is treated as an imposed stress coordinate, not a disease parameter.

See `docs/RAW_DATA_PHENOTYPES_V0_7_3.md` and `docs/RAW_DATA_PHENOTYPES_FREEZE_MANIFEST_2026-08-19.md`. Frozen code is under `pipelines/raw_data_phenotypes/`; compact results are under `results/raw_data_phenotypes/`.
