# Reproducibility tiers

## Tier 1: immediate figure reproduction

`code/reproduce_figures.sh` reads only committed CSV files and regenerates all
main figures plus the supporting robustness figures. It does not run stochastic
simulations.

## Tier 2: pipeline verification

`code/run_smoke_tests.sh` runs unit tests and reduced-data stochastic smoke
tests bundled with both pipelines. Smoke outputs are diagnostic only.

## Tier 3: full numerical reproduction

`code/run_full_analyses.sh` launches either full pipeline against
`data/inputs/`. The KL design comprises 264 scenarios, 32 dependent cell pairs,
31 path positions, five seeds, three integration steps, and two passes. The
nonequilibrium design comprises 264 scenarios, 31 path positions, and five
seeds at the primary integration step. Both pipelines support checkpoints.

Published outputs remain in `data/` and `results/validation/`; reruns write to
the corresponding pipeline directory and never overwrite the frozen release
tables automatically.
