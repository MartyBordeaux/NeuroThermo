# v1.3 outputs

Primary files:

- `component_boundary_summary_v1_3.csv` — scenario-first WT-exit, balance, and SCA3-entry curves for combined, kappa-only, and J-only surfaces.
- `biological_pair_component_boundaries_v1_3.csv` — within-pair uncertainty-aware boundary distributions and crossing-support weights.
- `component_effect_surface_summary_v1_3.csv` — biological-pair ensemble surface of kappa main effect, J main effect, combined effect, and non-additive interaction.
- `interaction_at_stage_boundaries_v1_3.csv` — component-effect magnitudes near the frozen ISI/active stage boundaries.
- `coupled_line_component_sensitivity_v1_3.csv` — local component sensitivities evaluated along the coupled line `q=p_intrinsic`.
- `RUN_SUMMARY.json` and `REPORT.md` — run audit and compact scientific summary.

Large/support files:

- `checkpoints/kappa_only/scenario_XXXX.csv.gz`
- `checkpoints/J_only/scenario_XXXX.csv.gz`
- `biological_pair_component_effect_surface_v1_3.csv.gz`
- `biological_pair_coupled_line_component_sensitivity_v1_3.csv`
- `biological_pair_interaction_at_stage_boundaries_v1_3.csv`
- `transition_pair_scenarios_v1_3.csv`

Figures:

- `01_boundary_WT_exit.png`
- `01_boundary_balance.png`
- `01_boundary_SCA3_entry.png`
- `02_interaction_surface_ISI.png`
- `03_coupled_line_component_sensitivity_ISI.png`
