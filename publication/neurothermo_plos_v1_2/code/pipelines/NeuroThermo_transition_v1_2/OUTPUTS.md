# Outputs

The full run writes `results_transition_ensemble_v1_2/`.

Key files:

- `ensemble_surface_summary.csv` — pair-balanced WT→SCA3 surface for all pairs and the core-secure subset.
- `biological_pair_surface_summary.csv.gz` — within-pair uncertainty-weighted surface for each of the 72 biological pairs.
- `transition_surface_states.csv.gz` — full support-scenario state table.
- `boundary_curve_summary.csv` — required drive at fixed intrinsic progress, and required intrinsic progress at fixed drive, for WT-exit/balance/SCA3-entry.
- `biological_pair_boundary_crossings.csv.gz` — pair-level boundary crossings before ensemble aggregation.
- `drive_sensitivity_surface.csv` — local derivatives of the primary `A_ISI` surface and the drive-dominance fraction.
- `drive_sensitivity_at_stage_boundaries.csv` — derivative summary near the three stage contours.
- `legacy_path_surface_recovery_summary.csv` — diagnostic comparison between the 2D surface and the frozen v1.1 early/coupled/late staging.
- `RUN_SUMMARY.json` — run counts and frozen-input hashes.

Figures:

- `01_primary_ISI_stage_map_core_secure.png`
- `02_required_drive_boundaries_core_secure.png`
- `03_drive_dominance_map_core_secure.png`
- `04_secondary_active_rate_stage_map_core_secure.png`
