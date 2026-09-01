# NeuroThermo

Frozen analysis repository for the manuscript *Multiscale dynamical staging in WT–SCA3 Purkinje-cell models*, covering Hindmarsh–Rose ensemble fitting, uncertainty-aware dynamical characterization, and WT→SCA3 transition staging of Purkinje-cell current-clamp recordings.

## Publication release

The self-contained publication bundle is in
[`publication/neurothermo_plos_v1_2/`](publication/neurothermo_plos_v1_2/README.md).
It contains analysis-ready data, complete KL-convergence and nonequilibrium
pipelines, Python and R figure code, publication figures, validation records,
and result-to-code traceability. All executable paths resolve within that
directory; no external server path is required.

## Current frozen stage

The active analysis is now:

1. **cell-fit v3.9** — final four-parameter HR cell fits;
2. **post-fit characterization v1.0** — animal-aware parameter and phenotype characterization;
3. **dynamic characterization v2.1** — experiment-support-restricted suprathreshold dynamics and near-optimal-solution robustness;
4. **endpoint ensemble v1.0.1** — uncertainty-aware WT and SCA3 endpoint clouds;
5. **transition ensemble v1.1** — corrected latency-invariant staging geometry;
6. **transition ensemble v1.2.1** — scenario-first intrinsic × drive map;
7. **transition ensemble v1.3.1** — factorized drive decomposition into model input scaling `kappa_I`, experimental current protocol `J`, and their interaction;
8. **KL convergence v1.0.1** — full-coverage multiseed convergence analysis;
9. **nonequilibrium geometry v1.0.1** — full-state Fisher geometry and detailed-balance audit.

### Frozen cell-fit model

- Fitted cell-level HR coordinates: `b`, `r`, `s`, `kappa_I`.
- Fixed HR constants: `a=1`, `c=1`, `d=5`, `x_R=-1.6`.
- Primary cohort: **18 multi-sweep cells (12 WT, 6 SCA3)**.
- Secondary descriptive set: 2 accepted single-sweep cells.
- Exact additive first-spike alignment per spiking sweep; no time rescaling and no last-spike anchoring.
- Binary rheobase bracket retained as the absolute excitability constraint.
- Final common bounds: `b=[0.5,7]`, `r=[1e-4,0.1]`, `s=[0.05,15]`, `kappa_I=[4e-4,4]`.

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

### Drive decomposition

The combined drive is decomposed into fitted input scaling `kappa_I`, applied
current density `J`, and their non-additive interaction. Along the coupled path,
the early combined-drive curve is J-shaped up to a parallel shift, whereas the
late curve is kappa-I-shaped. This supports a control handoff rather than a
single component dominating the entire path.

### Distributional and nonequilibrium results

Under full scenario, seed, grid, and integration-step coverage, full-state KL
balance occurs at a lower constructed path coordinate than firing-phenotype
balance at the ensemble level. This ordering is not interpreted as disease
time, causal precedence, or irreversible progression.

The slow coordinate contributes substantially to path Fisher geometry, and the
coarse time-reversal audit rejects detailed balance across the path. The
stationary object is therefore reported as a nonequilibrium potential, not as a
physical Hamiltonian.

## Repository layout

```text
publication/neurothermo_plos_v1_2/      self-contained publication release
docs/                                    frozen methods, scientific status, provenance
results/v3_9/                            final cell parameters
results/characterization_v1/             animal-aware characterization
results/dynamic_v2_1/                    support-restricted dynamic characterization
results/endpoint_ensemble_v1_0_1/        transition-ready endpoint ensemble
results/transition_v1_1/                 primary corrected staging
results/transition_v1_2_1/               scenario-first 2D intrinsic × drive analysis
results/transition_v1_3/                 drive decomposition and corrected figures
```

Earlier working pipelines remain recoverable from Git history. The publication
release is the authoritative entry point for reproducing article figures and
the two final stochastic analyses.
