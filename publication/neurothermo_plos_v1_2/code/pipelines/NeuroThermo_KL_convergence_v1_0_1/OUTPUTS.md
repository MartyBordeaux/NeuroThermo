# Outputs

The full run writes the following files to
`results_kl_convergence_v1_0_1/`.

- `KL_CONVERGENCE_VERDICT.json`: three-level manuscript decision, fatal versus
  supporting classification, and every gate.
- `RUN_SUMMARY.json`: design, task counts, fingerprint, and embedded verdict.
- `SCIENTIFIC_FINGERPRINT.txt`: checkpoint and result identity.
- `INPUT_MANIFEST.json`: hashes of the four frozen transition inputs.
- `scenario_markers.csv`: six marker variants for every scenario, step, and
  distributional view.
- `scenario_convergence_gates.csv`: time-step, seed, and method stability for
  every scenario and view.
- `pair_markers_both_orders.csv`: marker-first and curve-first pair estimates.
- `cell_balanced_deltas.csv`: endpoint-cell-balanced KL-minus-firing positions.
- `animal_pair_balanced_deltas.csv`: crossed animal-day summaries.
- `leave_one_animal_out.csv`: omission sensitivity for every animal/day.
- `ensemble_convergence_summary.csv`: manuscript-level summaries by step,
  view, marker, and aggregation order.
- `grid_retention.csv`: density-grid coverage diagnostics.
- `GRID_COVERAGE_AUDIT.json`: formal full-coverage verdict and task counts.
- `pilot_extents.csv`: per-task coordinate extrema from the pilot pass.
- `reference_grid_bounds.csv`: final bounds of every scenario-specific grid.
- `Fig_KL_convergence.pdf` and `.png`: convergence overview.
- `checkpoints/`: resumable path-level results; no raw trajectories.
- `extent_checkpoints/`: resumable pilot extrema; no raw trajectories.
- `reference_grids/`: fixed scenario-specific density grids.
