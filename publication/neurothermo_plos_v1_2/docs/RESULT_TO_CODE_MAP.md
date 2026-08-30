# Result-to-code map

| Article result | Source data | Executable code | Frozen check |
|---|---|---|---|
| Endpoint phenotype and model agreement | `data/figure_source/fig1_endpoint_cells.csv` | `code/figures/python/render_figures.py`, `code/figures/R/01_endpoint.R` | Figure 1 |
| Firing-phenotype staging along the constructed path | `data/figure_source/fig2_core_secure_curves.csv`, `fig2_primary_isi_staging.csv` | `code/figures/python/render_figures.py`, `code/figures/R/02_transition.R` | Figure 2 |
| Early current-dominated and late coupling-dominated combined drive | `data/figure_source/fig3_drive_surface_core_secure.csv`, `fig3_combined_drive_handoff_summary.csv` | `code/figures/python/render_figures.py`, `code/figures/R/03_decomposition.R` | Figure 3 |
| Full-state KL balance occurs at lower constructed p than firing balance | `data/kl_convergence_v1_0_1/` | `code/pipelines/NeuroThermo_KL_convergence_v1_0_1/`, `code/figures/python/render_fig4_multiseed.py` | `KL_CONVERGENCE_VERDICT.json`, Figure 4, Figure S3 |
| Slow-coordinate contribution and detailed-balance violation | `data/nonequilibrium_geometry_v1_0_1/` | `code/pipelines/NeuroThermo_nonequilibrium_geometry_v1_0_1/`, `code/figures/python/render_nonequilibrium_summary.py` | `FORMALISM_VERDICT.json`, Figure 5 |

The constructed coordinate `p` orders model states. It is not interpreted as
disease time, a causal trajectory, or evidence of irreversible one-way
progression.
