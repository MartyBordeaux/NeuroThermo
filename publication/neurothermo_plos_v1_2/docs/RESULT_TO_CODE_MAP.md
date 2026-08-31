# Result-to-code map

| Article result | Upstream generator chain | Figure-facing source data | Renderer | Frozen check / current gap |
|---|---|---|---|---|
| Endpoint phenotype and model agreement | `NeuroThermo_cell_fit_v3_9_frozen_exact` -> `NeuroThermo_characterization_v1_0` -> `NeuroThermo_dynamic_v2_1` -> `NeuroThermo_endpoint_ensemble_v1_0_1` | `data/figure_source/fig1_endpoint_cells.csv` and related Fig1 tables | `code/figures/python/render_figures.py`, `code/figures/R/01_endpoint.R` | Figure 1 source data are frozen; exact script that assembles the final Fig1 CSV from endpoint outputs still needs to be identified or reconstructed from frozen publication code |
| Firing-phenotype staging along the constructed path | endpoint ensemble -> `NeuroThermo_transition_v1_0` -> `NeuroThermo_transition_v1_1` | `data/figure_source/fig2_core_secure_curves.csv`, `fig2_primary_isi_staging.csv` | `code/figures/python/render_figures.py`, `code/figures/R/02_transition.R` | Figure 2 source data and frozen transition outputs are present; portable endpoint-to-v1.1 wiring remains to be validated |
| Early current-dominated and late coupling-dominated combined drive | transition v1.1 -> `NeuroThermo_transition_v1_2` -> `NeuroThermo_transition_v1_2_1` -> `NeuroThermo_transition_v1_3_1` | `data/figure_source/fig3_drive_surface_core_secure.csv`, `fig3_combined_drive_handoff_summary.csv` and related Fig3 tables | `code/figures/python/render_figures.py`, `code/figures/R/03_decomposition.R` | Figure 3 source data and compact frozen v1.2.1/v1.3 results are present; portable transition chain still needs clean-clone validation |
| Full-state KL balance occurs at lower constructed p than firing balance | final transition/support-state scenarios -> `NeuroThermo_KL_convergence_v1_0_1` | `data/kl_convergence_v1_0_1/` | `code/figures/python/render_fig4_multiseed.py` | `KL_CONVERGENCE_VERDICT.json`, Figure 4, Figure S3; pipeline is complete in release |
| Slow-coordinate Fisher contribution and detailed-balance violation | final transition/support-state scenarios -> `NeuroThermo_nonequilibrium_geometry_v1_0_1` | `data/nonequilibrium_geometry_v1_0_1/` | `code/figures/python/render_nonequilibrium_summary.py` | `FORMALISM_VERDICT.json`, Figure 5; pipeline is complete in release |

## Publication-critical upstream inputs

The canonical exact cell-fit executable is `code/pipelines/NeuroThermo_cell_fit_v3_9_frozen_exact/`. Its clean-clone input validation passes with 20 accepted cells, 18 primary multi-sweep cells, 113 spiking fit sweeps and 4884 selected spikes after overrides.

`NeuroThermo_characterization_v1_0` requires the complete v3.9 `cell_fit_summary.csv` and `sweep_fit_summary.csv`; the current compact `results/upstream_frozen/v3_9/primary_cell_parameters.csv` is not sufficient to rerun characterization.

`NeuroThermo_dynamic_v2_1` requires the exact six-file frozen layer:

- `primary_cell_master.csv`
- `accepted_spiking_sweeps.csv`
- `selected_spike_events.csv`
- `threshold_brackets.csv`
- `final_identifiability_alternatives.csv`
- `animal_id_map.csv`

The publication animal-ID source table is now stored at `data/animal_id_recovery/accepted_cohort.csv`; a publication wrapper may convert it to the historical input format without changing scientific content.

The constructed coordinate `p` orders model states. It is not interpreted as disease time, a causal trajectory, or evidence of irreversible one-way progression.
