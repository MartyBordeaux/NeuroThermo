# Frozen design

The scientific design is frozen before the full server calculation.

- Cohort: `core_secure_all_support`, exactly 264 scenarios and 32 pairs.
- Path: coupled interpolation at 31 positions.
- Steps: 0.05, 0.025, 0.0125 ms.
- Seeds: 20260818, 21260821, 22260823, 23260837, 24260855.
- Primary step: 0.025 ms.
- Density: 22 bins per dimension. For each scenario, a pilot pass measures the
  observed minima and maxima over all three steps and all five seeds. The
  common grid spans those extrema plus a 2% margin, with a one-bin Gaussian
  filter and pseudocount 1e-10. Every analysis path must retain 100% of its
  samples.
- Views: `xyz`, `xy`, and `z`. Fatal gates protect the primary `xyz`
  seed-median isotonic result; `xy`, alternative markers, and animal-day
  summaries are supporting sensitivity analyses.
- Crossing persistence: three path positions.
- Numerical thresholds are unchanged from v1.0.0. Fatal versus supporting
  gate roles are declared in v1.0.1 before the full-coverage rerun; this
  hierarchy was introduced after diagnosing the overly conjunctive v1.0.0
  decision rule.
- No individual support scenario is treated as a biological replicate.
