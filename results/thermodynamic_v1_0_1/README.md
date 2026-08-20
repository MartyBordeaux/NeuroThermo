# Thermodynamic transition v1.0.1 — compact frozen results

This directory freezes the compact scientific outputs of the complete 264/264 support-scenario thermodynamic run and its prespecified diffusion/timestep robustness checks. Bulky per-scenario `.npz` checkpoints are intentionally excluded from GitHub.

Primary interpretation:

- KL information geometry is the primary robust stochastic-state progression result.
- Fisher information is a broad, heterogeneous sensitivity/reorganization measure, not a unique tipping point.
- Entropy is descriptive and does not define a robust transition stage.
- Model EPR supports nonequilibrium restructuring, but the exact EPR-peak location is diffusion-dependent.
- Exact finite-state Hatano–Sasa identity passes to numerical precision.
- Empirical Hatano–Sasa sampling remains rare-event limited.
- Crooks/path-IFT empirical gates do not pass; classical Jarzynski/work-Crooks remain disabled.

Readable working source is under `code/thermodynamic_transition_v1_0_1/`.
