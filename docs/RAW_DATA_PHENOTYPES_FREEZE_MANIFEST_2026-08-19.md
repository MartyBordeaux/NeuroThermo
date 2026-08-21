# Raw-data phenotype freeze manifest — 2026-08-19

## Frozen code snapshots

| repository path | role | status |
|---|---|---|
| `pipelines/raw_data_phenotypes/v0_3_1/` | base ABF-to-phenotype extraction | frozen |
| `pipelines/raw_data_phenotypes/v0_5_1/` | dependency-reduced endpoint fingerprint | frozen |
| `pipelines/raw_data_phenotypes/v0_6_1/` | calibrated current-resolved vulnerability | frozen |
| `pipelines/raw_data_phenotypes/v0_7_0/` | adjusted predictive-dynamics validation | frozen |
| `pipelines/raw_data_phenotypes/v0_7_1/` | adjacent-window shuffle-surrogate analysis | superseded but retained |
| `pipelines/raw_data_phenotypes/v0_7_2/` | non-overlap shuffle analysis and IAAFT attempt | shuffle result frozen; IAAFT excluded |
| `pipelines/raw_data_phenotypes/v0_7_3/` | rank-Gaussian exact-spectrum Fourier validation | frozen primary validation |

## Compact results

Selected aggregate exact tests, fidelity diagnostics and run summaries are stored under `results/raw_data_phenotypes/`. Native analysis manifests are omitted because they contain absolute server paths and private provenance hashes. Record- and cell-level tables, source paths, large surrogate-long tables, ABF recordings, generated figures, virtual environments and intermediate calibration exports are also intentionally omitted. The frozen implementation includes Python source, tests, dependency metadata and synthetic-data generation code. Per-version server configs, launch scripts and internal method notes are omitted because they expose the private server data layout; the repository-level report provides the sanitized method record.

## Superseded stages

Versions v0.2.x and v0.3.0 were calibration and extraction-development stages superseded by v0.3.1. v0.5.0 and v0.6.0 were corrected by v0.5.1 and v0.6.1. The v0.4.x animal-aware branch is not part of the frozen inference because the available animal count was inadequate and the project explicitly returned to cell-level inference.

## Reproducibility boundary

Re-running the code snapshots requires the original ABF files and private upstream result locations, which are not committed to GitHub. v0.7.2 has no completed final manifest because its IAAFT fidelity gate failed; that failure is documented rather than repaired retrospectively.
