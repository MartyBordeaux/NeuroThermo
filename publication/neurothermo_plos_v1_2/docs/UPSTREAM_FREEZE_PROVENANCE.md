# Upstream freeze provenance

Audit date: 2026-08-31

This file records the immutable project artifacts restored into the publication snapshot and the evidence used to identify the final versions. It distinguishes original frozen source from later portability wrappers.

## Readable pipeline source restored from `main`

| Publication path | Original project path | Git tree SHA | Role |
|---|---|---|---|
| `code/pipelines/NeuroThermo_cell_fit_v3_9/` | `code/cell_fit_v3_9/` | `4e51c49a9b19dd149ffd42c7fce61509d5284c57` | final restricted HR wide-bound cell refit and identifiability |
| `code/pipelines/NeuroThermo_characterization_v1_0/` | `code/characterization_v1_0/` | `6dc0d92266657cd6b22bb3ea051a6901daa12312` | post-fit characterization |
| `code/pipelines/NeuroThermo_transition_v1_3_1/` | `code/transition_v1_3_1/` | `f88188fecf525cead41abfdde48ca42a25f46db3` | factorized kappa_I/J decomposition and v1.3.1 plotting correction |

These trees were copied by Git object identity, not reconstructed or edited. Machine-specific configurations contained inside them therefore remain historical frozen source. Publication-portable configurations must be added separately.

## Compact frozen upstream results restored from `main`

| Publication path | Git tree SHA |
|---|---|
| `results/upstream_frozen/v3_9/` | `5d29f5b9f1cad2866fe3672213fef06aeedcbb27` |
| `results/upstream_frozen/characterization_v1/` | `464417790203b823aa1ef6b0e7db06cbd2193265` |
| `results/upstream_frozen/dynamic_v2_1/` | `af62bbe88929155217bcd5b35e46261a691b3748` |
| `results/upstream_frozen/endpoint_ensemble_v1_0_1/` | `d5cf53d92c4fcc98b847476c64730cace8bb01cc` |
| `results/upstream_frozen/transition_v1_1/` | `daa28db70336a12fe63ff165ee46ca78de3a5c47` |
| `results/upstream_frozen/transition_v1_2_1/` | `8a52788cb5e462d0d6c7f82c3eebc315189fb761` |
| `results/upstream_frozen/transition_v1_3/` | `805745944bcf3d365ccbe2c08bb05dd49df7daeb` |

These are compact publication-critical tables and summaries, not substitutes for the executable pipelines that generated them.

## Frozen package hashes recorded before transition analysis

From `docs/FREEZE_MANIFEST_2026-08-16.md` in the private analysis repository:

| Artifact | SHA-256 |
|---|---|
| `neurothermo_dynamic_v2.1.0.zip` | `f4fa72124600616b8246ec6d8ecd73982e0d65e4f0aa6f0fa720141bb1ffcae7` |
| `results_dynamic_v2_1(1).zip` | `1c942191a7c6b00ef12938e3d584b156afec5fecd283d2168afc3aaee4819fbc` |
| `neurothermo_endpoint_ensemble_v1.0.1.zip` | `448f0290a196bda9d57bc67eb3850e7b2778ecc7bdce83eae7e1939b828d9c87` |
| `results_endpoint_ensemble_v1_0(2).zip` | `17e05a8955a19bca69a7453a1b6bde02aa3946e6acada182efe3770528e12aec` |

## Frozen transition package/result hashes

From `docs/TRANSITION_FREEZE_MANIFEST_2026-08-18.md`:

| Artifact | SHA-256 |
|---|---|
| `results_transition_ensemble_v1_1(1).zip` | `e0116bf9a3e935f1e79f4cffe44219d1f3cb63641485183e7ec7d620c6962453` |
| `results_transition_ensemble_v1_2_1(1).zip` | `a96b54a1a3ceb3602db2dbbd24368f486b2d31e6b6873db601a01bcea502bc5b` |
| `results_transition_ensemble_v1_3.zip` | `d64028132a921c3adb0916d2eab7df6083e24ec18245650cac1fe75ac2dd5745` |
| `results_transition_ensemble_v1_3_corrected.zip` | `02951b60eed7709de2753379c38d1238455293db79fa89fafb366efa473ecd97` |
| `neurothermo_transition_ensemble_v1.3.1.zip` | `74b0ef218306c3543248e8e34f3c408abf1ead7320c7cbb73f126995d4fda8ff` |

The corrected v1.3 result archive differs from the original numerical archive only by the three corrected boundary figures and `FIGURE_FIX.md`; v1.3 numerical results are unchanged.

## Unresolved calibration-input provenance

The final checked-in `cell_fit_v3_9` configuration requires:

- `candidate_events_with_predictions.csv` with `fixed_qc_detected` inclusion;
- `frozen_accepted_spiking_sweeps_v3_6.csv`;
- `frozen_peak_overrides_v3_6.csv`;
- `frozen_threshold_brackets_v3_6.csv`;
- v3.1 baseline cell/sweep/identifiability tables;
- `seed_cell_summary_v3_9.csv`.

The current publication requirement had separately listed v3.5-named frozen files. Exact searches of the accessible GitHub trees and File Library have not yet located either complete v3.5 or v3.6 calibration bundles. Therefore the publication snapshot must not claim that v3.5 and v3.6 files are identical and must not rename either version. This remains the principal blocker to an end-to-end clean-clone rerun beginning before v3.9 fitting.

## Encoded historical executable packages

Branch `freeze-working-code-through-v1.3.1` preserves encoded package artifacts for `dynamic_v2_1`, `endpoint_v1_0_1`, `transition_v1_1`, `transition_v1_2_1`, and related stages. Their existence establishes recoverability, but the publication release requires decoded, readable source trees before it can be considered complete. The SHA-256 package hashes above are the acceptance references for decoded recovery.
