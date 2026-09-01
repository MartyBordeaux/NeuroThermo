# Outputs

`results_transition_ensemble_v1_1/` contains:

- `PRIMARY_ISI_STAGING.csv` — primary stage locations by path family in the core-secure biological-pair subset.
- `transition_projection_reference_v1_1.csv` — frozen ISI reference plus corrected experimental active-rate reference.
- `transition_projection_transform_v1_1.csv` — robust scaling parameters.
- `active_rate_reference_cells.csv` — the secure endpoint cells that define the corrected active-rate boundary.
- `transition_paths_reprojected.csv` — the original 121,524 simulated states with `active_support_rate_hz`, corrected `A_active`, and independently recomputed `A_isi_v1_1`.
- `scenario_stage_markers_v1_1.csv` — scenario-level ISI, active-rate, and dual-agreement markers.
- `biological_pair_stage_summary_v1_1.csv` — parameter-support uncertainty collapsed within each biological WT×SCA3 pair.
- `path_family_stage_summary_v1_1.csv` — all-pair and core-secure-pair stage summaries.
- `biological_pair_curve_summary_v1_1.csv`, `ensemble_curve_summary_v1_1.csv` — pair-balanced curves.
- `occupancy_at_primary_ISI_stages.csv` — occupancy at the ISI-defined WT-exit, balance, and SCA3-entry positions.
- `CORRECTION_AUDIT.csv` — numerical checks that no simulated value was changed and that ISI reprojection is identical to v1.0.
- `RUN_SUMMARY.json`.
- `figures/` — compact primary staging and projection-agreement figures.
