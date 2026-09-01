# v3.9.0

- Built directly on v3.8.
- Changed only the lower bound of `s`: 0.25 -> 0.05 for all cells.
- Retained v3.8 bounds for `b`, `r`, `kappa_I`, exact first-spike alignment, rheobase constraint, objective, QC, and identifiability separation criterion.
- Refits the full accepted cohort under one common parameter space.
- Uses final v3.8 solutions as seeds/baselines.
- Added `s_boundary_stress_summary.csv` and explicit lower-bound occupancy diagnostics.
- No special branch or sensitivity analysis for SCA3_05.
