# Computational DAG audit — publication release

Audit date: 2026-08-31

Target branch: `publication/plos-v1.2`

Scientific reference for this audit: the final manuscript and Supporting Information dated 2026-08-29, plus the frozen analysis documentation on `main`.

## Release acceptance criterion

A clean clone of the eventual public repository must be sufficient to install the documented environment and regenerate the reported numerical tables, source-data tables, and figures without editing paths or supplying files from `/root/neurothermo`, `~/neurothermo`, another Git branch, or an unpublished server directory.

The release must preserve the frozen restricted Hindmarsh–Rose formulation: fitted `b`, `r`, `s`, `kappa_I`; fixed `a=1`, `c=1`, `d=5`, `x_R=-1.6`; exact additive first-spike alignment; no time rescaling; no last-spike anchoring; nonspiking sweeps used only through the binary rheobase bracket; `SCA3_05` retained in the accepted fit cohort.

## Current release defect

`code/run_full_analyses.sh` currently reruns only `kl` and `nonequilibrium`. Figures 1–3 and S1 are rendered from already prepared CSV source tables, so the release cannot currently regenerate those tables from the frozen electrophysiology/fit inputs. The release is therefore a figure/source-data snapshot plus two complete stochastic pipelines, not yet a complete publication workflow.

## Required computational DAG

Status values:

- `COMPLETE_RELEASE`: executable code and required inputs are already present in the publication release.
- `RECOVERABLE`: final code/results exist elsewhere in the private repository or frozen Git history and must be restored into the release.
- `NONPORTABLE`: code exists but has machine-specific or unpublished-file dependencies that must be removed.
- `MISSING`: required artifact/source has not yet been located in the accessible GitHub snapshot.
- `SCOPE_CHECK`: requested analysis exists in project history but is not referenced by the 2026-08-29 main/SI manuscript; retain separately until article linkage is established.

| ID | Stage | Final source/version | Direct publication role | Current status | Required remediation |
|---|---|---|---|---|---|
| 00 | Raw current-clamp + provenance/QC | archived ABF + QC/provenance workbook | provenance and reproducibility of spike/QC inputs | MISSING | add shareable/de-identified source data or an explicit external-data accession plus machine-readable provenance/QC table |
| 01 | Event/QC calibration and frozen spike selections | frozen accepted-sweep manifest, manual peak overrides, threshold brackets, event table | defines the exact fit target | MISSING / provenance conflict | recover exact frozen files and hashes; resolve `v3_5` versus `v3_6` naming/content from provenance rather than renaming files |
| 02 | Restricted HR cell fit | `cell_fit_v3_9` | cell-specific endpoint parameters, fit loss, spike-fit diagnostics, identifiability layer | RECOVERABLE + NONPORTABLE | copy readable package into release; replace `/root/neurothermo/...` and `../calibration/...`; bundle all required calibration inputs and seeds; preserve full v3.9 outputs, not only `primary_cell_parameters.csv` |
| 03 | Post-fit characterization | `characterization_v1_0` | endpoint parameter/phenotype summaries and parameter associations | RECOVERABLE | restore code; replace manual ZIP/XLSX arguments with release-relative paths; preserve generated machine-readable tables |
| 04 | Experiment-support-restricted dynamics | `dynamic_v2_1` | Fig. 1/S1 firing-rate and mean-ISI support tables | RECOVERABLE | decode/restore frozen package from `freeze-working-code-through-v1.3.1`; wire inputs to v3.9/characterization outputs |
| 05 | Endpoint uncertainty ensemble | `endpoint_ensemble_v1_0_1` | uncertainty-aware endpoint support states and core-secure definition | RECOVERABLE | decode/restore frozen package; preserve cell weights/support-state weights and endpoint source tables |
| 06 | Corrected transition staging | `transition_v1_1` | Fig. 2 path-family WT-exit/balance/SCA3-entry locations | RECOVERABLE | decode/restore final v1.1 package and inputs; do not substitute earlier transition versions |
| 07 | Intrinsic × combined-drive surface | `transition_v1_2_1` | Fig. 3A/B and scenario-first uncertainty aggregation | RECOVERABLE | decode/restore v1.2.1 package and full scenario inputs/results |
| 08 | Factorized input decomposition | `transition_v1_3_1` | Fig. 3C/D: combined, `kappa_I`, applied `J`, interaction | RECOVERABLE | copy readable package from `main`; make v1.2 input path release-relative; preserve v1.3 numerical results and v1.3.1 figure correction |
| 09 | Earlier stochastic thermodynamic overlay | `thermodynamic_transition_v1_0_1` | supporting provenance for Hatano–Sasa/rare-event development; not source of final full-coverage KL verdict | RECOVERABLE | archive as supporting pipeline if retained; clearly separate from final KL/nonequilibrium branches |
| 10 | Full-coverage endpoint-relative KL | `NeuroThermo_KL_convergence_v1_0_1` | main Fig. 4; S2/S3; Tables S5–S7 | COMPLETE_RELEASE | retain exact package, configs, seeds, validation fingerprints and full-coverage outputs |
| 11 | Nonequilibrium geometry + Fisher/time-reversal audit | `NeuroThermo_nonequilibrium_geometry_v1_0_1` | slow-coordinate FI result; detailed-balance audit; supporting tables | COMPLETE_RELEASE | retain exact package/configs/validation; keep continuous-current EPR excluded after failed stationarity gate |
| 12 | Publication source-data assembly | current `data/figure_source/*` | all figure-facing source tables | PARTIAL | add scripts that generate every source table from upstream stage outputs; precomputed CSVs may remain as frozen checks but not as the only source |
| 13 | Figure/table rendering | Python/R scripts in `code/figures` | main and supporting figures/tables | PARTIAL | retain renderers; add explicit upstream source-table generation and output checksum validation |
| 14 | Master workflow | current `run_full_analyses.sh` | reviewer entry point | PARTIAL | replace two-target dispatcher with full ordered workflow plus stage-specific commands, resume rules, and smoke profile |
| 15 | Environment lock | mixed `requirements.txt` files | software reproducibility | PARTIAL | consolidate exact Python/R versions; add lock/freeze output compatible with the documented Python runtime |
| 16 | Pathwise temporal-order analysis | project-requested analysis | requested publication archive component | MISSING / SCOPE_CHECK | locate final package/results; map to a specific current manuscript/SI claim before inserting into the article workflow |
| 17 | 2D KL auxiliary analysis | project-requested analysis | requested supporting analysis | MISSING / SCOPE_CHECK | locate final code/results and distinguish it from the full-state `xyz` and fast `xy` marginals already implemented in KL v1.0.1 |
| 18 | PI/Fourier analysis | project-requested analysis | requested publication archive component | MISSING / SCOPE_CHECK | locate final package/results; current 2026-08-29 HR manuscript does not use predictive-information/Fourier results, so keep provenance separate until linkage is explicit |

## Confirmed upstream chain for the 2026-08-29 HR manuscript

```text
raw/QC/event selections
    -> cell_fit_v3_9
    -> characterization_v1_0
    -> dynamic_v2_1
    -> endpoint_ensemble_v1_0_1
    -> transition_v1_1
    -> transition_v1_2_1
    -> transition_v1_3_1
    -> publication source tables for endpoint/transition figures

endpoint/transition support-state scenarios
    -> KL_convergence_v1_0_1
    -> Fig. 4 + KL supporting figures/tables

endpoint/transition support-state scenarios
    -> nonequilibrium_geometry_v1_0_1
    -> Fisher/time-reversal supporting results
```

`thermodynamic_transition_v1_0_1` is retained as a separate earlier stochastic overlay unless a direct final-output dependency is demonstrated. Its Hatano–Sasa and rare-event diagnostics must not be silently treated as the source of the final full-coverage KL or Fisher/time-reversal results.

## Blocking provenance issue: frozen v3.5 versus v3.6 inputs

The current publication requirement names:

- `frozen_accepted_spiking_sweeps_v3_5.csv`
- `frozen_peak_overrides_v3_5.csv`
- `frozen_threshold_brackets_v3_5.csv`

The checked-in final `cell_fit_v3_9` server configuration instead requests `frozen_accepted_spiking_sweeps_v3_6.csv`, `frozen_peak_overrides_v3_6.csv`, and `frozen_threshold_brackets_v3_6.csv`, plus several v3.1 baseline/seed files. This discrepancy is publication-blocking until file contents, hashes, and provenance are recovered. No file will be renamed merely to satisfy the expected version label.

## Immediate remediation order

1. Recover and hash the exact frozen QC/calibration inputs used by v3.9.
2. Restore readable final packages for dynamic v2.1, endpoint v1.0.1, transition v1.1, and transition v1.2.1 from the frozen Git branch.
3. Copy cell-fit v3.9, characterization v1.0, transition v1.3.1, and any retained thermodynamic-transition package into the publication tree.
4. Replace all machine-specific paths with release-root-relative configuration paths.
5. Add source-table assembly scripts so Figures 1–3/S1 no longer depend only on precomputed CSVs.
6. Replace `run_full_analyses.sh` with a complete staged dispatcher and add a single clean-clone `run_all` route.
7. Freeze Python and R environments and record random seeds/configuration hashes.
8. Validate outputs against the frozen publication CSVs/figures using checksums and scientific-value checks.
9. Only after the private publication snapshot passes clean-clone validation, copy the snapshot to a separate public repository.
