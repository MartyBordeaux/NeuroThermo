# NeuroThermo

Frozen analysis repository for Hindmarsh–Rose ensemble fitting and uncertainty-aware dynamical characterization of Purkinje-cell current-clamp recordings in WT and SCA3.

## Current frozen stage

The active analysis is now:

1. **cell-fit v3.9** — final four-parameter HR cell fits;
2. **post-fit characterization v1.0** — animal-aware parameter and phenotype characterization;
3. **dynamic characterization v2.1** — experiment-support-restricted suprathreshold dynamics and near-optimal-solution robustness;
4. **endpoint ensemble v1.0.1** — uncertainty-aware WT and SCA3 endpoint clouds prepared for the next WT→SCA3 transition/staging analysis.

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

The primary robust phenotype is higher capacitance-normalized rheobase in SCA3, lower firing rate in SCA3 at experiment-supported `q=0.75`, and longer mean ISI in SCA3 at the same supported state.

### Transition-ready endpoint representation

Endpoint ensemble v1.0.1 contains **64 actual HR support states** from 18 biological cells. Near-optimal solutions are retained as within-cell model uncertainty and do not receive independent biological weight.

The primary transition projection is

`(log10 J_rheo, log10 firing_rate_q75)`

with mean ISI at `q=0.75` retained as a parallel validation coordinate rather than discarded as redundant.

## Repository layout

```text
docs/                              frozen methods, scientific status, provenance
data/animal_id_recovery/            recovered animal assignments
results/v3_9/                       final cell parameters
results/characterization_v1/        animal-aware post-fit characterization
results/dynamic_v2_1/               support-restricted dynamic characterization
results/endpoint_ensemble_v1_0_1/   uncertainty-aware transition-ready endpoints
```

Earlier working pipelines remain recoverable from Git history. The current branch is intentionally a compact frozen-analysis snapshot rather than a chronological dump of intermediate pipelines.

See `docs/SCIENTIFIC_STATUS.md`, `docs/DYNAMIC_V2_1.md`, `docs/ENDPOINT_ENSEMBLE_V1_0_1.md`, `docs/METHODS_V3_9.md`, and `docs/ANIMAL_ID_PROVENANCE.md`.
