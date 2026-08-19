# Changelog

## 0.7.2

- Replaced adjacent ordinal-code pairs with a prespecified non-overlapping 4 ms primary lag.
- Added 8, 16 and 32 ms lag sensitivity.
- Replaced classic AAFT with iterative AAFT.
- Added per-surrogate spectrum, log-PSD, autocorrelation, amplitude and convergence diagnostics.
- Added frozen IAAFT fidelity gates.
- Renamed `PI_excess` to the signed `surrogate_centered_PI`.
- Added lower-tail, upper-tail and two-sided surrogate p-values.
- Added exact maxT correction across lags and across lag×current combinations.
- Retained all 13 WT and 7 SCA3 cells with cell-level inference only.
