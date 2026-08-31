# v2.1 outputs

- `observed_current_dynamics_all_solutions.csv` — best and near-optimal models evaluated at the same actually recorded J; includes experimental descriptors and support-restricted model descriptors.
- `observed_current_phase_profiles_best.csv` — best-fit phase profiles at actual J, restricted to cycles inside the experimental first-to-last-spike interval.
- `observed_current_scalar_robustness.csv` — best-vs-alternative scalar robustness at identical J.
- `observed_current_phase_robustness.csv` — phase-profile robustness at identical J.
- `rheobase_refinement_all_solutions.csv` — independent rheobase refinement for best and alternative parameter sets.
- `q_support_by_cell.csv` — observed q range of every primary cell and support flags for requested q targets.
- `q_interpolated_experiment.csv` — experimental spike descriptors interpolated only within observed-current support.
- `q_interpolated_model_all_solutions.csv` — model descriptors at the same supported q targets; no extrapolation.
- `q_scalar_robustness_near_optimal.csv` — robustness of q-interpolated scalar phenotype.
- `q_phase_profiles_best.csv` — best-fit supported phase profiles at q targets.
- `q_phase_robustness_near_optimal.csv` — robustness of phase profiles after support-restricted q interpolation.
- `group_q_medians.csv` — cell-level descriptive group medians, with experiment and best model kept separate.
- `animal_q_medians.csv` — descriptive medians by recovered animal ID.
- `group_q_phase_median_profiles.csv` — group phase profiles from best fits at supported q targets.
- `RUN_SUMMARY.json` — frozen counts, q coverage, and compact robustness diagnostics.

The 20% robustness thresholds are diagnostic flags only; they are not inferential significance thresholds.
