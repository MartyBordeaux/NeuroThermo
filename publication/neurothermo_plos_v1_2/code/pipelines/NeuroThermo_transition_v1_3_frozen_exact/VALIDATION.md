# Validation

The package is designed around the frozen transition v1.2 result archive/directory and requires all 988 v1.2 combined-surface checkpoints.

Validation checks include:

1. source `RUN_SUMMARY.json` reports transition v1.2.0;
2. all 988 combined checkpoints exist;
3. the full server profile selects exactly 988 support-state scenarios;
4. WT/SCA3 biological-pair weights are inherited unchanged from the frozen scenario table;
5. at `p_component = p_intrinsic`, kappa-only and J-only states reduce to the same coupled parameter/current state as the combined map, up to numerical rheobase tolerance;
6. J-only rheobase is computed once per intrinsic row because changing experimental J does not change the HR parameter vector or its model rheobase.

The primary ISI stage boundaries are frozen at:

- WT-exit: `A_ISI = 0.1358293233470019`;
- balance: `A_ISI = 0.5`;
- SCA3-entry: `A_ISI = 0.7978563373093712`.

The analysis uses scenario-first uncertainty propagation for boundary crossings, matching the corrected v1.2.1 logic.
