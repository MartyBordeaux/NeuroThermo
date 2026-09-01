# Frozen design

- Model and diffusion: identical to the multiseed KL ensemble v1.0.0.
- Cohort: `core_secure_all_support`.
- Expected inputs: 264 supported parameter scenarios representing 32 dependent
  WT/SCA3 cell-pair combinations.
- Path: 31 equally spaced `p` values.
- Seeds: 20260818, 21260821, 22260823, 23260837, 24260855.
- Integrator: Euler–Maruyama, `dt=0.025 ms`.
- Per stationary run: 2400 ms burn-in, 6000 ms sampling, saved every 0.5 ms.
- Diffusion: `(0.0025, 0.01, 0.00025)`.
- Primary density: 22 bins per state dimension, Gaussian smoothing 1 bin.
- Markov coarse-graining: fixed 4 × 4 × 3 quantile partition, 5 ms primary lag.
- Fluctuation schedules: 15 strictly increasing indices selected from the same
  31-state path for linear, xyz-FI, xy-FI, and friction protocols.
- Animal sensitivity: four recoverable WT animals and two recoverable SCA3
  animal-days in the selected core-secure cells.

Nothing in this package modifies the endpoint cohort, interpolation rule, noise,
or stationary simulation settings used by the frozen multiseed analysis.

All inputs use the final units reported in the manuscript.
