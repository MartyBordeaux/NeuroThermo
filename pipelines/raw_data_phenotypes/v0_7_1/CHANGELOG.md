# Changelog

## 0.7.1

- Recomputes ordinal predictive information from raw ABF stationary intervals.
- Adds a hard reproduction gate against the frozen v0.3.1 PI values.
- Adds 500 value-shuffled surrogates per sweep as the primary finite-sample/overlap null.
- Adds 500 AAFT surrogates per sweep as a secondary spectral sensitivity analysis.
- Freezes shuffle-median-subtracted PI cell-AUC as the primary endpoint.
- Enumerates all 77,520 cell-label assignments.
- Adds current-wise maxT tests and surrogate-count convergence.
- Repeats label-blind activity/technical residualization and secondary nested WT burden scoring with shuffle-corrected PI.
- Retains all 20 cells and removes animal-level inference.
- Keeps `I_exit` secondary and explicitly excludes disease-time interpretation.
