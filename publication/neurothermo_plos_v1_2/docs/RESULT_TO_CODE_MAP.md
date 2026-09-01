# Result-to-code map

| Article result | Upstream generator chain | Figure-facing source data | Canonical renderer | Clean-clone status |
|---|---|---|---|---|
| Endpoint phenotype and model agreement | raw ABF -> `NeuroThermo_stage1_qc_fixed` -> frozen accepted sweeps/overrides -> `NeuroThermo_cell_fit_v3_9_frozen_exact` -> `NeuroThermo_characterization_v1_0` -> `NeuroThermo_dynamic_v2_1` -> `NeuroThermo_endpoint_ensemble_v1_0_1` -> `code/assemble_figure_source.py` | `data/figure_source/fig1_endpoint_cells.csv` | `code/figures/python/render_figures.py` | PASS. Fig1 source reconstruction matches the frozen table to <=5.6e-17. Capacitance-aware publication transformation is explicit in the assembler. Unresolved WT recording labels are derived from raw ABF header dates and independently validated. |
| Firing-phenotype staging along the constructed path | endpoint ensemble -> `NeuroThermo_transition_v1_0` -> `NeuroThermo_transition_v1_1` -> `code/assemble_figure_source.py` | `fig2_core_secure_curves.csv`, `fig2_primary_isi_staging.csv`, `fig2_projection_reference.csv` | `code/figures/python/render_figures.py` | PASS. v1.1 reprojection and v1.1->v1.2 assembly pass; Fig2 source reconstruction matches frozen source to <=1.5e-14. |
| Early current-dominated and late coupling-dominated combined drive | transition v1.1 -> `NeuroThermo_transition_v1_2` -> `NeuroThermo_transition_v1_2_1` -> `NeuroThermo_transition_v1_3_frozen_exact` -> `code/assemble_figure_source.py` | `fig3_drive_surface_core_secure.csv`, `fig3_coupled_component_sensitivity_core_secure.csv`, `fig3_drive_sensitivity_at_boundaries_core_secure.csv`, `fig3_interaction_at_boundaries_core_secure.csv` | `code/figures/python/render_figures.py` | PASS. Portable v1.2/v1.2.1/v1.3 validation passes across all 988 v1.2 checkpoints; four renderer-consumed Fig3 tables reconstruct to <=4.5e-16. `fig3_combined_drive_handoff_summary.csv` is auxiliary and is not consumed by the canonical Fig3 renderer. |
| Full-state KL balance occurs at lower constructed p than firing balance | final transition/support-state scenarios -> `NeuroThermo_KL_convergence_v1_0_1` | `data/kl_convergence_v1_0_1/` | `code/figures/python/render_fig4_multiseed.py` | COMPLETE. `KL_CONVERGENCE_VERDICT.json`, Figure 4 and supporting convergence outputs are included. |
| Slow-coordinate Fisher contribution and detailed-balance violation | final transition/support-state scenarios -> `NeuroThermo_nonequilibrium_geometry_v1_0_1` | `data/nonequilibrium_geometry_v1_0_1/` | `code/figures/python/render_nonequilibrium_summary.py` | COMPLETE. `FORMALISM_VERDICT.json` and the nonequilibrium/Fisher summary figure are included. This is the publication Fisher-information implementation; no separate Fisher package is required by the current manuscript. |

## Publication-facing source-data assembly

`code/assemble_figure_source.py` is the canonical Fig1--3 source-data assembler. It writes recomputed tables to `results/recomputed/figure_source/` and compares them against the immutable release references in `data/figure_source/`; it never overwrites those references.

The assembler currently reconstructs all eight Fig1--3 tables consumed by the canonical Python renderer:

- `fig1_endpoint_cells.csv`
- `fig2_core_secure_curves.csv`
- `fig2_primary_isi_staging.csv`
- `fig2_projection_reference.csv`
- `fig3_drive_surface_core_secure.csv`
- `fig3_coupled_component_sensitivity_core_secure.csv`
- `fig3_drive_sensitivity_at_boundaries_core_secure.csv`
- `fig3_interaction_at_boundaries_core_secure.csv`

For Fig1, the publication capacitance-aware transformation is explicit: capacitance and `kappa_I` are doubled, rheobase current density `J=I/Cm` is halved, and rheobase in pA is unchanged. Recovered animal-day labels come from the canonical animal map; the four unresolved WT recordings use the date in their raw ABF header as a recording-provenance label only, not as a recovered animal identity.

For Fig2, the core-secure cohort is rebuilt from endpoint `core_q75_secure`: 8 WT x 4 SCA3 endpoint cells -> 32 dependent biological pairs, which are then aggregated with equal pair weighting. The same capacitance-aware correction halves the displayed rheobase current density.

## Publication-critical upstream inputs

The canonical exact cell-fit executable is `code/pipelines/NeuroThermo_cell_fit_v3_9_frozen_exact/`. Its clean-clone input validation passes with 20 accepted cells, 18 primary multi-sweep cells, 113 spiking fit sweeps and 4884 selected spikes after overrides.

`NeuroThermo_characterization_v1_0` uses the complete imported v3.9 result tree in `data/v3_9_results_full/`; the compatibility XLSX required by the historical reader is generated deterministically from `data/animal_id_recovery/accepted_cohort.csv`.

`NeuroThermo_dynamic_v2_1` uses the exact six-file frozen layer: `primary_cell_master.csv`, `accepted_spiking_sweeps.csv`, `selected_spike_events.csv`, `threshold_brackets.csv`, `final_identifiability_alternatives.csv`, and `animal_id_map.csv`.

The constructed coordinate `p` orders model states. It is not interpreted as disease time, a causal trajectory, or evidence of irreversible one-way progression.

## Renderer policy

The canonical reproducible figure renderer is Python/matplotlib under the pinned Python 3.9.25 environment. The R scripts are retained as optional alternative plotting code and are not required for reviewer clean-clone reproduction.
