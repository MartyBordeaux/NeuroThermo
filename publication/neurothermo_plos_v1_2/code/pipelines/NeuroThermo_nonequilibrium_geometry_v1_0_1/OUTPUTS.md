# Output map

| File | Meaning |
|---|---|
| `FORMALISM_VERDICT.json` | NESS/equilibrium-candidate decision and Jarzynski–Crooks gate |
| `PROTOCOL_VERDICT.json` | Whether adaptive schedules improve both error and ESS at equal state count |
| `NUMERICAL_VALIDATION.json` | FI, KL, chain-rule, current-QC, and schedule-uniqueness checks |
| `RUN_SUMMARY.json` | Completion counts, input location, fingerprint, preflight status |
| `INPUT_MANIFEST.json` | SHA-256 hashes of the four frozen input tables |
| `preflight_fixed_points.csv` | Equilibria, Jacobian spectra, local stability |
| `preflight_current_continuation.csv` | Deterministic current continuation and oscillation classification |
| `preflight_endpoint_membership.csv` | Every scenario assigned to its WT and SCA3 endpoint groups |
| `stationary_geometry.csv` | Scenario × seed × p primary estimators |
| `pair_balanced_geometry.csv` | Support-weighted scenario results within each dependent cell pair |
| `ensemble_geometry.csv` | Cell-pair-balanced medians and IQRs |
| `animal_pair_balanced_geometry.csv` | Equal-cell-pair means within animal pairs |
| `animal_balanced_geometry.csv` | Equal-animal-pair medians and IQRs; primary robustness summary |
| `animal_mapping_used.csv` | Auditable cell-to-animal mapping for the selected cohort |
| `local_kl_fisher_check.csv` | Adjacent KL versus Fisher quadratic approximation |
| `markov_cycle_affinities.csv` | Supported three-state cycle affinities |
| `thermodynamic_length_summary.csv` | Full path lengths by scenario and seed |
| `adaptive_protocols.csv` | Constant-speed schedules based on xyz-FI, xy-FI, and friction |
| `pair_balanced_protocols.csv` | Support-balanced schedules within cell pairs |
| `animal_pair_balanced_protocols.csv` | Schedules aggregated within animal pairs |
| `potential_marginals.csv` | Relative one-dimensional potential marginals at selected p |
| `fluctuation_relations.csv` | Exact/sampled Hatano–Sasa and sampled path-ratio IFT diagnostics |
| `pair_balanced_fluctuation_relations.csv` | Support-balanced fluctuation diagnostics |
| `animal_pair_balanced_fluctuation_relations.csv` | Animal-pair sensitivity for fluctuation diagnostics |
| `protocol_performance_summary.csv` | Fair schedule comparison with equal unique-state counts |
| `oscillatory_endpoint_pair_balanced_geometry.csv` | Sensitivity restricted to scenarios with oscillatory endpoints |
| `markov_cache/` | Compressed transition matrices and invariant distributions for cheap protocol reanalysis |
| `estimator_sensitivity.csv` | Density-grid and Markov-lag sensitivity analyses |
| `failures.csv` | Failed scenario-seed tasks, if continuation after error is enabled |
| `figures/` | Three PNG/PDF diagnostic figures |

Checkpoint JSON and Markov cache files are scientific-fingerprint protected. A
config or frozen-input change invalidates stale checkpoints automatically.
