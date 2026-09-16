# Upstream path map for the revision pipelines

The revision pipeline configs are preserved exactly as run and therefore contain original server paths. On the publication branch, the corresponding frozen upstream data are inherited from `publication/neurothermo_plos_v1_2/`.

Core path mapping:

| Original server input | Publication-package counterpart |
|---|---|
| `.../neurothermo_transition_ensemble_v1_2/results_transition_ensemble_v1_2` | `../neurothermo_plos_v1_2/data/transition_v1_2_results/` |
| `.../neurothermo_endpoint_ensemble_v1_0/results_endpoint_ensemble_v1_0` | `../neurothermo_plos_v1_2/data/endpoint_ensemble_v1_0_results/` |
| `.../neurothermo_transition_ensemble_v1_3/results_transition_ensemble_v1_3` | `../neurothermo_plos_v1_2/data/transition_v1_3_results/` |
| v3.5 frozen fit/QC inputs | v1.2 calibration bundle / deterministic calibration extraction |
| raw WT/SCA3 ABFs | `../neurothermo_plos_v1_2/data/raw/WT/` and `../neurothermo_plos_v1_2/data/raw/SCA3/` |
| recovered animal/cell mapping | `../neurothermo_plos_v1_2/data/animal_id_recovery/accepted_cohort.csv` |

The v1.2 package README documents the clean-clone commands that verify and extract the calibration bundle. Run those preflight/calibration steps before attempting a full revision-pipeline recomputation outside the original server layout.

For the chain of new analyses:

1. run `neurothermo_route_robustness_v1_0_0` from the v1.2 endpoint/transition result trees;
2. point `neurothermo_current_intervention_predictions_v1_0_3` to the generated route results plus the v1.2 frozen QC/raw inputs;
3. point `neurothermo_raw_proxy_robustness_v1_0_0` to the generated v1.0.3 prediction results plus the same frozen raw/QC inputs.

The exact server configs are retained for provenance instead of silently rewriting the paths after analysis.
