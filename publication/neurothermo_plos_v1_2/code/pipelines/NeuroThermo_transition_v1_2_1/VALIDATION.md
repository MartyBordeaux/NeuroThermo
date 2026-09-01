# Validation logic

Validation requires the frozen v1.2 result set: 988 checkpoint surfaces, 72 biological pairs, 32 core-secure pairs and 31×31 grids. The sum of scenario weights must equal one.

The correction does not recompute HR trajectories. A crossing is first detected on each scenario surface. Within each biological pair, marker q25/median/q75 are then computed with the same finite-marker weighted-quantile convention used in v1.1. Parameter-support states with no marker are not silently discarded: their contribution is reported separately as `crossing_support_weight`.

The old one-dimensional drive-early/coupled/drive-late paths are reconstructed by bilinear interpolation of every scenario's 2D surface and compared with frozen v1.1 ISI staging.

## Reference full post-processing run

The supplied v1.2 result archive was processed without any new HR simulations:

- 988 support-state scenarios;
- 72 biological WT×SCA3 pairs;
- 32 core-secure pairs;
- 31×31 intrinsic×drive grid inherited from v1.2;
- 949,468 frozen state rows represented by the checkpoint set.

A clean reference run takes about 30–35 s in the build environment. Drive-early is recovered from the 2D surface essentially exactly. Coupled and late SCA3-entry are close to frozen v1.1. The larger residual error for drive-late WT-exit is retained explicitly as a 31×31 grid-resolution diagnostic and is not tuned away.
