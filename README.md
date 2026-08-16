# NeuroThermo

Frozen analysis repository for Hindmarsh–Rose ensemble fitting of Purkinje-cell current-clamp recordings in WT and SCA3.

## Current frozen stage

The active analysis is **cell-fit v3.9 + post-fit characterization v1.0**.

- Four fitted cell-level HR coordinates: `b`, `r`, `s`, `kappa_I`.
- Fixed HR constants: `a=1`, `c=1`, `d=5`, `x_R=-1.6`.
- Primary cohort: **18 multi-sweep cells (12 WT, 6 SCA3)**.
- Secondary descriptive set: 2 accepted single-sweep cells.
- Exact additive first-spike alignment per spiking sweep; no time rescaling and no last-spike anchoring.
- Binary rheobase bracket is retained as the absolute excitability constraint.
- Final common bounds: `b=[0.5,7]`, `r=[1e-4,0.1]`, `s=[0.05,15]`, `kappa_I=[2e-4,2]`.
- Identifiability is evaluated in the same wide parameter space with the separation criterion inherited from the narrower reference range rather than inflated with the wider bounds.

## Repository layout

```text
docs/                         frozen methodology and scientific interpretation
data/animal_id_recovery/       recovered animal assignments and provenance
results/v3_9/                  final cell parameters and identifiability
results/characterization_v1/   animal-aware descriptive characterization
```

The previous per-cell/alternative-model working tree has been removed from the current branch to prevent confusion with the frozen analysis. It remains available through Git history.

## Scientific status

The most stable SCA3-associated HR coordinate is increased `s`, with a more moderate increase in `b`. The apparent WT–SCA3 difference in `r` does not survive wide-bound stress and `r` is the least identifiable coordinate. Raw `kappa_I` strongly covaries with capacitance and is not interpreted as an independent SCA3 phenotype. Capacitance-normalized rheobase is markedly higher in the recorded SCA3 cells.

Recovered animal identifiers show only **two recovered animals per group** in the current primary dataset. Cells are therefore not treated as independent biological replicates for genotype-level inference. Cell-level effects are descriptive for the recorded-cell ensemble; animal-level summaries are reported as medians without formal p-values.

See `docs/SCIENTIFIC_STATUS.md`, `docs/METHODS_V3_9.md`, and `docs/ANIMAL_ID_PROVENANCE.md`.
