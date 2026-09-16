# Revised result-to-code map

| Revised manuscript element | Primary code | Committed source/audit data | Upstream dependency |
|---|---|---|---|
| Experimental endpoints and HR fit adequacy | historical v1.2 pipelines | `endpoint_plot_data.csv`, `archival_endpoint_scalar_summary.csv` | `../neurothermo_plos_v1_2/` |
| 2D full-state KL vs firing-balance geometry | `neurothermo_route_robustness_v1_0_0` | `contour_ordering_by_pair.csv`, `contour_ordering_summary.csv`, `discordance_region_summary.csv` | v1.2 endpoint + transition result trees |
| 525-route robustness, including 512 random monotone routes | `neurothermo_route_robustness_v1_0_0` | `route_family_summary.csv`, `ROUTE_ROBUSTNESS_VERDICT.json` | same |
| Pair/scenario robustness and aggregation checks | `neurothermo_route_robustness_v1_0_0` | `selected_scenarios.csv`, `pairwise_interactions_summary.csv`, `ROUTE_ROBUSTNESS_VERDICT.json` | same |
| Factorial/Shapley supporting attribution | `neurothermo_route_robustness_v1_0_0` | `shapley_summary.csv` | same |
| Direct physical-current predictions at fixed SCA endpoint parameters | `neurothermo_current_intervention_predictions_v1_0_3` | `direct_current_predictions_by_cell.csv`, `FROZEN_P1_P2_TARGETS.csv`, `CURRENT_PREDICTION_VERDICT.json` | route results + v1.2 raw/QC/frozen inputs |
| Retrospective current-prediction checks | `neurothermo_current_intervention_predictions_v1_0_3` | `retrospective_prediction_tests.csv`, `RAW_VALIDATION_REPORT.json`, `v1_0_3_primary_reproduction.csv` | same archival cells |
| Raw voltage-proxy sensitivity/robustness audit | `neurothermo_raw_proxy_robustness_v1_0_0` | `cell_robustness_classification.csv`, `prediction_robustness_summary.csv`, `reference_jackknife_summary.csv`, `bootstrap_summary.csv`, `RAW_PROXY_ROBUSTNESS_VERDICT.json` | current-prediction results + v1.2 raw/QC inputs |
| Main text and Supporting Information | manuscript source artifact | all files above | v1.2 + revision layer |

Notes:

- Cross-combinations, support scenarios, stochastic seeds, grid nodes, and routes are computational units, not biological replicates.
- The joint model-space drive coordinate is not a direct physical-current intervention. The current-only pipeline is the separate intervention calculation.
- Protocol current `J` in factorial attribution is an experimenter-controlled factor and is not interpreted as a disease mechanism.
- Negative raw-proxy findings are intentionally retained because they delimit the biological interpretation of the model-space result.
