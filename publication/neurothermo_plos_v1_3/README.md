# NeuroThermo manuscript revision package v1.3

This directory aligns the repository with the 14 September 2026 manuscript draft:

**Dynamical state reorganization precedes the firing phenotype in WT–SCA3 Purkinje-cell model space**

It is a revision layer on top of the frozen reproducibility release `v1.2` (commit `245791ac47af49feb9b765407831c24d0cdf6b60`). The v1.2 package remains the authoritative source for raw ABF recordings, QC/provenance inputs, frozen spike selections, the four-parameter Hindmarsh–Rose fit, endpoint ensembles, and the pre-existing transition/KL analyses.

## What this revision layer records

The revised draft additionally depends on three post-rejection analyses:

- `neurothermo_route_robustness_v1_0_0`: 2D intrinsic–drive atlas, 525-route robustness audit, numerical convergence checks, and factorial attribution;
- `neurothermo_current_intervention_predictions_v1_0_3`: fixed-SCA-endpoint current-only predictions and retrospective raw-waveform checks;
- `neurothermo_raw_proxy_robustness_v1_0_0`: frozen sensitivity audit of the retrospective voltage proxy.

The exact source-archive and result-archive SHA-256 fingerprints are recorded in `CHECKSUMS.sha256` and `docs/FULL_RESULT_ARCHIVES.md`. `docs/RESULT_TO_CODE_MAP.md` maps the revised manuscript claims to those analyses and their final source tables. `docs/UPSTREAM_PATH_MAP.md` maps the original server paths in the frozen configs to the upstream v1.2 publication package.

## Scientific scope

The revised manuscript's principal result is a model-space statement: within the prescribed endpoint-anchored intrinsic–drive surface, full-state KL affinity and the reduced firing phenotype define separated balance boundaries. The separation is robust across the frozen route ensemble, including 512 random monotone routes. The route coordinate is not disease time and route count is not biological replication.

The direct-current and raw-proxy analyses are prediction/audit layers. Their mixed or negative retrospective support is retained in the final interpretation; they are not presented as independent validation of the model-space ordering.

## Deliberate exclusions

Post-draft exploratory branches that did not yield a manuscript-level result are not part of this manuscript-aligned layer: Hatano–Sasa route calculations, RQA/Takens analyses, direct-relaxation/capacitance-staging analyses, the 50-pA-only audit, and the state–relaxation bridge. None is cited or used by the 14 September draft.

## Archival status

The historical v1.2 package is fully present in this repository. For the revised analyses, this branch currently records their exact identities, provenance, result-to-code map, and cryptographic fingerprints. The complete post-rejection source/result ZIP archives must still be attached to the permanent archival release before the manuscript's final Data and Code Availability placeholder can be replaced by a permanent release identifier/DOI.

This distinction is intentional: the repository must not claim that a binary archive is committed when only its verified fingerprint is present.

## Historical metadata not inferred

The original animal-ethics approval identifier/coverage and the sex distribution of the recovered cohort still require external historical confirmation. They are not inferred computationally.
