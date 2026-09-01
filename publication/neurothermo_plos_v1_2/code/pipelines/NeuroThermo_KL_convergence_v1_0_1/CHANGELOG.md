# Changelog

## v1.0.1

- Replaced the one-seed quantile grid with a two-pass all-task-extrema grid.
- Enforced 100% retained sample mass for every analysis path.
- Added resumable pilot-extent checkpoints and formal coverage outputs.
- Replaced the all-or-none verdict with prespecified fatal and supporting gate
  tiers while retaining every original numerical threshold.
- Added full-coverage and verdict-tier tests.

The stochastic model, frozen cohort, path, time steps, seeds, stationary-run
lengths, density bin count, KL definition, marker variants, aggregation
procedures, and numerical gate thresholds are unchanged from v1.0.0.
