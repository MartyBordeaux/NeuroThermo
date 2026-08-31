# Validation rules

1. Frozen cohort must contain 18 primary cells: 12 WT and 6 SCA3.
2. Endpoint support must contain 64 admissible HR states and enumerate 988 WT×SCA3 support scenarios across 72 biological pairs.
3. Scenario weights must sum to one before any profile filtering.
4. q=.75 support must exist for all 18 primary cells.
5. The primary v1.1 ISI boundaries must reproduce `0.135829` and `0.797856` to 1e-4.
6. The 2D interpolation must reproduce exact WT at `(0,0)` and exact SCA3 at `(1,1)` for every support scenario.
7. Changing only `p_drive` must change only `kappa_I` and `J`; changing only `p_intrinsic` must change `b,r,s` and the support window while leaving `kappa_I,J` fixed.
8. The three legacy early/coupled/late curves are diagnostics embedded in the surface, not independent biological hypotheses.
9. `drive` in v1.2 remains the combined `(kappa_I,J)` protocol coordinate. The planned `kappa_I`/`J` decomposition is a later stage and must not be inferred from v1.2 alone.

## Local execution checks

The package was tested with:

- `pytest`: 4 tests passed;
- smoke profile: 2 support scenarios × 7×7 grid;
- preliminary profile: all 72 best×best biological pairs × 21×21 grid = 31,752 simulated states;
- server-profile validation: 988 support scenarios × 31×31 grid = 949,468 planned states.

The preliminary run generated all summaries and four figures successfully. The full all-support server profile was not run locally because it is intentionally the VPS workload.
