# Computational DAG audit — publication release

Audit date: 2026-08-31

Target branch: `publication/plos-v1.2`

Scientific reference: final manuscript and Supporting Information dated 2026-08-29 plus frozen analysis documentation.

## Acceptance criterion

A clean clone of the eventual public repository must install the documented environment and regenerate the reported numerical tables, source-data tables and figures without editing paths or supplying files from `/root/neurothermo`, `~/neurothermo`, another Git branch or an unpublished server directory.

The frozen restricted Hindmarsh–Rose formulation is preserved: fitted `b`, `r`, `s`, `kappa_I`; fixed `a=1`, `c=1`, `d=5`, `x_R=-1.6`; exact additive first-spike alignment; no time rescaling; no last-spike anchoring; nonspiking sweeps only through the binary rheobase bracket; `SCA3_05` retained in the accepted cohort.

## Major provenance findings resolved on 2026-08-31

1. The exact server calibration archive was imported into the publication branch with SHA-256 `0c930506021826aec8ee2987fe83cd4a1537fa42b6d3fad335a5520fcbb610bd`.
2. The v3.5- and v3.6-named accepted-sweep, peak-override and threshold-bracket files are byte-identical. The publication layer therefore uses the canonical v3.5 names while preserving the original archive and hashes.
3. The historical exact frozen `cell_fit_v3_9` package was recovered from immutable Git blobs. It differs from the later readable `main` copy in 14 of 24 files, including core loader/objective/optimizer/pipeline code. The exact frozen package is now the canonical publication executable.
4. A clean GitHub-hosted runner using Python 3.9.25 successfully reconstructed and hash-checked all calibration files and validated the exact v3.9 input layer: 20 accepted cells, 18 primary multi-sweep cells, 113 spiking fit sweeps, 4884 spikes after overrides, 8 peak overrides and 20 threshold brackets.
5. Readable exact/frozen source is now present for `dynamic_v2_1`, `endpoint_ensemble_v1_0_1`, transition v1.0, v1.1, v1.2, v1.2.1, and exact frozen v3.9. `transition_v1_3_1`, characterization, KL and nonequilibrium source are also present.

## Current blocker

The release can now start the exact v3.9 analysis reproducibly, but the downstream chain is not yet fully executable from the publication tree because full v3.9 outputs and the prepared frozen input layer used by `dynamic_v2_1` are not yet bundled. `characterization_v1_0` requires full `cell_fit_summary.csv` and `sweep_fit_summary.csv`, not only `primary_cell_parameters.csv`. `dynamic_v2_1` requires a six-file frozen layer: `primary_cell_master.csv`, `accepted_spiking_sweeps.csv`, `selected_spike_events.csv`, `threshold_brackets.csv`, `final_identifiability_alternatives.csv`, and `animal_id_map.csv`.

A self-hosted runner audit confirmed that the runner account `github-runner` has no passwordless sudo and cannot inspect `/root/neurothermo` directly. Server artifacts therefore need a controlled export into a runner-readable location or must be regenerated from the publication workflow.

## Stage status

| ID | Stage | Final source/version | Status | Publication action |
|---|---|---|---|---|
| 00 | Raw current-clamp + QC/provenance | archived ABF + QC workbook | MISSING_PUBLICATION_LAYER | add redistributable/de-identified raw/QC package or external accession |
| 01 | Frozen event/QC calibration | exact server archive + canonical v3.5 names | COMPLETE_RELEASE | retained with SHA-256 verification |
| 02 | Restricted HR fit | `NeuroThermo_cell_fit_v3_9_frozen_exact` | INPUT_VALIDATED; FULL_OUTPUTS_NOT_BUNDLED | full rerun is now launchable; preserve/reproduce complete v3.9 outputs |
| 03 | Characterization | `NeuroThermo_characterization_v1_0` | SOURCE_PRESENT / INPUT_BLOCKED | provide full v3.9 outputs and publication animal map |
| 04 | Support-restricted dynamics | `NeuroThermo_dynamic_v2_1` | SOURCE_PRESENT / FROZEN_INPUT_LAYER_MISSING | recover or deterministically generate six-file frozen layer |
| 05 | Endpoint uncertainty ensemble | `NeuroThermo_endpoint_ensemble_v1_0_1` | SOURCE_PRESENT / INPUT_WIRING_PENDING | wire to dynamic results and validate |
| 06 | Transition base | `NeuroThermo_transition_v1_0` | SOURCE_PRESENT / INPUT_WIRING_PENDING | wire endpoint support states |
| 07 | Corrected transition staging | `NeuroThermo_transition_v1_1` | SOURCE_PRESENT / INPUT_WIRING_PENDING | replace historical relative-path discovery with publication wrapper |
| 08 | Intrinsic × drive transition | `NeuroThermo_transition_v1_2` + `v1_2_1` | SOURCE_PRESENT / INPUT_WIRING_PENDING | wire v1.1→v1.2→v1.2.1 explicitly |
| 09 | Factorized drive | `NeuroThermo_transition_v1_3_1` | SOURCE_PRESENT / INPUT_WIRING_PENDING | point to publication v1.2 results |
| 10 | Full-coverage KL | `NeuroThermo_KL_convergence_v1_0_1` | COMPLETE_RELEASE | retain exact configs/seeds/results |
| 11 | Nonequilibrium geometry/Fisher | `NeuroThermo_nonequilibrium_geometry_v1_0_1` | COMPLETE_RELEASE | retain exact configs/seeds/results |
| 12 | Figure source-data assembly | current `data/figure_source/*` | PARTIAL | map/generate Fig1–3 source tables from upstream outputs |
| 13 | Figure rendering | Python/R renderers | PRESENT | retain and add source checksum validation |
| 14 | Master workflow | `code/run_full_analyses.sh` | PARTIAL | currently supports calibration/preflight/cellfit validation plus KL/nonequilibrium; add downstream stages |
| 15 | Environment lock | Python 3.9.25 + mixed requirements | PARTIAL | freeze exact transitive Python and R versions |
| 16 | Pathwise temporal order | requested project component | SCOPE_CHECK | locate and map to current manuscript/SI claim before article DAG inclusion |
| 17 | 2D KL auxiliary | requested project component | SCOPE_CHECK | distinguish from existing `xy` marginal in KL v1.0.1 |
| 18 | PI/Fourier | requested project component | SCOPE_CHECK | locate; current final HR manuscript linkage not established |

Machine-readable details are in `PUBLICATION_WORKFLOW_INVENTORY.tsv`.

## Confirmed computational chain

```text
raw/QC
  -> exact frozen calibration
  -> exact cell_fit_v3_9
  -> characterization_v1_0
  -> dynamic_v2_1 frozen input assembly
  -> dynamic_v2_1
  -> endpoint_ensemble_v1_0_1
  -> transition_v1_0
  -> transition_v1_1
  -> transition_v1_2
  -> transition_v1_2_1
  -> transition_v1_3_1
  -> Fig1–3/S1 source-data assembly

transition/support-state scenarios
  -> KL_convergence_v1_0_1
  -> Fig4 + KL supplementary tables/figures

transition/support-state scenarios
  -> nonequilibrium_geometry_v1_0_1
  -> Fisher/time-reversal supplementary results
```

`thermodynamic_transition_v1_0_1` remains a separate earlier stochastic overlay unless a direct dependency of the final publication outputs is demonstrated.

## Exact calibration hashes

- archive: `0c930506021826aec8ee2987fe83cd4a1537fa42b6d3fad335a5520fcbb610bd`
- candidate events: `af35c327b313482f534aa59669a47e52a4078f912a5e342efcfddf0158455640`
- accepted sweeps v3.5: `dad46b831eb4613af4a49673f83854e4ef48b81d0934c087234562d81a447a54`
- peak overrides v3.5: `64e35808199e6108355b015b4ca9ded6070deed852927877e705ccf118e95069`
- threshold brackets v3.5: `47ba271e6b8d70704de1c49aaac3677c6ee21e3001f33faaafad8761177f9741`
- v3.1 cell summary: `85b1fa2c457e4affc0db438cf885b4406f61b943cbc08073fbdeb7f4b57f42f9`
- v3.1 sweep summary: `5663d59c35aeb105ee45b0c4c8606375210294f377a6ee3adcd771356a70ab12`
- v3.1 identifiability: `16e810e3331a0f6eb6bc1c815bb0e0d5574ee93966b95c00346522d5470957d1`
- v3.9 seed summary: `cb74bc0783c9fd1db11cacba13ccabd273cfc225e6dc019ab6e4215433dceb72`

## Next remediation order

1. Recover or regenerate the full v3.9 output directory needed by characterization and dynamic input assembly.
2. Recover the exact dynamic six-file frozen input layer or document/implement its deterministic generation from v3.9 + animal-map artifacts.
3. Add publication wrappers for characterization, dynamic, endpoint and transition stages; leave frozen scientific source unchanged.
4. Validate each stage against the frozen result tables already committed under `results/upstream_frozen/`.
5. Map Fig1–3 source tables to their generating stage outputs and update `RESULT_TO_CODE_MAP.md`.
6. Freeze exact Python/R dependency versions and seeds.
7. Resolve raw/QC redistribution and project-requested scope-check analyses.
8. Only after a complete clean-clone run passes, create the separate public repository snapshot.
