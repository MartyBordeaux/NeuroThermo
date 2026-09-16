# NeuroThermo manuscript revision package v1.3

This directory aligns the public repository with the 14 September 2026 manuscript draft:

**Dynamical state reorganization precedes the firing phenotype in WT–SCA3 Purkinje-cell model space**

It is a revision layer on top of the frozen reproducibility release `v1.2` (commit `245791ac47af49feb9b765407831c24d0cdf6b60`). The v1.2 package remains the source for raw ABF recordings, QC/provenance inputs, frozen spike selections, the four-parameter Hindmarsh–Rose fit, endpoint ensembles, and the pre-existing transition/KL analyses. This v1.3 directory adds only the analyses and manuscript material required by the revised draft.

## Added for the revised draft

- `artifacts/neurothermo_route_robustness_v1_0_0.zip`: 2D intrinsic–drive atlas, 525-route robustness audit, numerical convergence checks, and factorial attribution.
- `artifacts/neurothermo_current_intervention_predictions_v1_0_3.zip`: fixed-SCA-endpoint current-only predictions and retrospective raw-waveform checks.
- `artifacts/neurothermo_raw_proxy_robustness_v1_0_0.zip`: frozen sensitivity audit of the retrospective voltage proxy.
- `artifacts/neurothermo_final_publication_source_data_compact.zip`: compact publication-facing verdicts and tables used to write and audit the revised Results and Supporting Information.
- `artifacts/NeuroThermo_manuscript_source_2026-09-14.zip`: source of the 14 September 2026 draft and Supporting Information.
- `docs/AUDIT_2026-09-16.md`: repository-to-manuscript audit.
- `docs/RESULT_TO_CODE_MAP.md`: map from revised manuscript claims to code and source data.
- `docs/FULL_RESULT_ARCHIVES.md`: checksums for the complete result archives produced by the three revision pipelines.

## Scientific scope

The revised manuscript's principal result is a model-space statement: within the prescribed endpoint-anchored intrinsic–drive surface, full-state KL affinity and the reduced firing phenotype define separated balance boundaries. The separation is robust across the frozen route ensemble, including 512 random monotone routes. The route coordinate is not disease time and route count is not biological replication.

The direct-current and raw-proxy analyses are retained as prediction/audit layers. Their mixed or negative retrospective support is part of the final interpretation; they are not presented as independent validation of the model-space ordering.

## What is intentionally not included in the manuscript-aligned revision layer

Post-draft exploratory branches that did not yield a manuscript-level result are not part of this release layer, including the later Hatano–Sasa route experiment, RQA/Takens exploration, direct-relaxation/capacitance-staging experiments, and the state–relaxation bridge. Their exclusion is deliberate rather than a missing dependency: none is cited or used by the 14 September draft.

## Reproduction relationship to v1.2

Clone the `publication/plos-v1.3` branch. The full historical upstream data are already present in `publication/neurothermo_plos_v1_2/` on the same branch. The revision pipeline archives retain the exact server configs used for the reported calculations; original absolute server paths are preserved for provenance. `docs/UPSTREAM_PATH_MAP.md` identifies their corresponding locations in the v1.2 publication package.

The complete large result archives are listed by SHA-256 in `docs/FULL_RESULT_ARCHIVES.md`. Compact manuscript-facing summaries are included in the committed source-data artifact so the reported values can be audited without the large archives.

## Known metadata items not invented by this repository update

The manuscript still requires external confirmation of the original animal-ethics approval identifier/coverage and the sex distribution of the recovered cohort. These are historical experimental metadata and are not inferred computationally.
