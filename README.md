# NeuroThermo

Reproducibility repository for the manuscript revision **Dynamical state reorganization precedes the firing phenotype in WT–SCA3 Purkinje-cell model space**.

## Manuscript-aligned publication material

The historical clean-clone reproducibility release remains:

- `publication/neurothermo_plos_v1_2/`
- tag `v1.2`
- frozen commit `245791ac47af49feb9b765407831c24d0cdf6b60`

It contains the raw ABF recordings, QC/provenance inputs, frozen manual selections, the restricted four-parameter Hindmarsh–Rose fits, endpoint/transition analyses, expensive frozen result trees, and clean-clone replay infrastructure.

The September 2026 manuscript revision adds a separate alignment layer:

- `publication/neurothermo_plos_v1_3/`

That directory documents the post-rejection analyses used by the revised manuscript: the 2D/525-route robustness analysis, direct-current prediction analysis, and raw voltage-proxy robustness audit. It also records the exact pipeline/result archive checksums and a manuscript result-to-code map. The v1.2 package remains the authoritative upstream source for raw data and frozen historical analyses.

## Current scientific scope

The primary result is a conditional model-space statement. Within the prescribed endpoint-anchored intrinsic–drive Hindmarsh–Rose surface, full-state KL affinity and the reduced rheobase–ISI firing phenotype define separated balance boundaries. The separation persists across a frozen ensemble of 525 routes, including 512 random monotone routes. Route coordinates are not interpreted as disease time, and route/scenario counts are not treated as biological replication.

The final retrospective physical-current and raw-voltage analyses are retained as prediction/audit layers. Their mixed or negative support limits biological interpretation rather than validating the KL boundary experimentally.

Post-draft exploratory branches that did not produce a manuscript-level result—Hatano–Sasa, RQA/Takens, relaxation/capacitance staging, the 50-pA relaxation audit, and the state–relaxation bridge—are intentionally excluded from the manuscript-aligned release layer.

## Core frozen model/data facts

- Fitted HR parameters: `b`, `r`, `s`, `kappa_I`.
- Fixed HR constants: `a=1`, `c=1`, `d=5`, `x_R=-1.6`.
- Primary multi-sweep cohort: 18 cells (12 WT, 6 SCA3).
- Core-secure model-space subset: 8 WT + 4 SCA3 cells.
- Dependent WT×SCA3 endpoint combinations: 32.
- Retained support scenarios: 264.
- Frozen route ensemble: 525 routes, including 512 random monotone routes.
- SCA3 recordings map to two recovered animal-day groups; no population-level genotype inference is claimed.

See `publication/neurothermo_plos_v1_3/docs/AUDIT_2026-09-16.md` for the manuscript-to-repository audit and `publication/neurothermo_plos_v1_3/docs/FULL_RESULT_ARCHIVES.md` for exact archive fingerprints.
