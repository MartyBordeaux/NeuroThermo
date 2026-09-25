# Hatano-Sasa: production core-32 analysis

This folder contains the current production Hatano-Sasa analysis and replaces the earlier WT_04–SCA3_04 pilot files.

## Scope

- 8 core-secure WT cells × 4 core-secure SCA3 cells = 32 dependent WT–SCA3 computational contrasts.
- Best-fit layer only. Retained near-optimal parameter variants are not treated as independent biological replicates and are not propagated in this script.
- Exact q = 0.75 endpoint current is reconstructed from the frozen cell-specific rheobase and observed spiking-current support:
  J(q) = Jrheo + q (Jmax_observed - Jrheo).
- Coupled path: b, s, and J linear; r and kappa_I logarithmic.
- 31 path positions, stochastic dt = 0.025 ms, D = (0.0025, 0.01, 0.00025).
- 12 s burn-in and 6 s stationary sampling.
- Training and independent holdout seed families are disjoint.
- Hatano-Sasa is evaluated in a 48-state (4×4×3) finite-state coarse graining.
- The script also checks deterministic x(t) spiking along every path and summarizes local Y around the coarse KL-balance region.

## Main result files

- `production_pair_summary_primary.csv` — one row per WT×SCA3 pair with the primary independent-start HS results and QC metrics.
- `sca3_stratified_HS_summary.csv` — HS performance stratified by SCA3 endpoint.
- `Y_KL_region_comparison.csv` — local Y behavior in early, KL-window, pair-KL-window, and late regions.
- `WT04_SCA304_corrected_q75_snapshot.csv` — corrected q=0.75 reference result for the original WT_04–SCA3_04 pair.
- `core_endpoints_best_q75.csv` — frozen endpoint parameters and exact observed-current support used to reconstruct q=0.75.
- `hs_core32_discrete_crossfit_production.py` — reproducible production script.

## Reproduce

Install:

```bash
python -m pip install numpy pandas scipy numba
```

Run from this folder:

```bash
python hs_core32_discrete_crossfit_production.py --workers 8
```

The default output directory is `results_core32_discrete_q75/`.

## Interpretation boundary

The exact finite-state Hatano-Sasa identity is a mathematical consistency check of the estimated stationary Markov representation. The scientifically relevant control is the independent-start result. Across the 32 pairs, deviations are concentrated in endpoints with extremely slow z dynamics (notably SCA3_06 and SCA3_09), so stationarity/mixing limitations must remain explicit in interpretation.
