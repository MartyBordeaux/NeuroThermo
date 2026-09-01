# Computational DAG audit — publication release

Audit date: 2026-08-31

Target branch: `publication/plos-v1.2`

Scientific reference: current PLOS Computational Biology manuscript result map plus frozen analysis documentation.

## Acceptance criterion

A clean clone of the eventual public repository must install the documented environment and regenerate or validate the reported numerical tables, source-data tables and figures without editing paths or supplying files from `/root/neurothermo`, `~/neurothermo`, another Git branch or an unpublished server directory.

The frozen restricted Hindmarsh–Rose formulation is preserved: fitted `b`, `r`, `s`, `kappa_I`; fixed `a=1`, `c=1`, `d=5`, `x_R=-1.6`; exact additive first-spike alignment; no time rescaling; no last-spike anchoring; nonspiking sweeps only through the binary rheobase bracket; `SCA3_05` retained in the accepted cohort.

## Current validated state

The earlier provenance and portability blockers have been resolved in the publication tree.

1. The release contains all 50 raw current-clamp ABF recordings plus the QC/provenance workbook and SHA-256 raw-data manifest.
2. The exact server calibration archive is embedded and hash-verified. Canonical frozen v3.5 accepted-sweep, peak-override and threshold-bracket files are materialized reproducibly.
3. Clean-clone raw-to-QC recomputation on the pinned historical environment reproduces 6800 candidate events, 6217 classifier-positive events, 6039 final fixed-QC spikes and 186 manually changed events, with zero missing/extra identities and zero decision mismatches. Thirteen classifier probabilities differ only at machine-precision scale (maximum absolute difference `1.1102230246251565e-16`) and do not cross the decision threshold.
4. The exact frozen `cell_fit_v3_9` source is the canonical publication executable. The complete v3.9 result tree required by downstream stages is embedded via the verified upstream server bundle and materialized by `code/prepare_upstream_inputs.py`.
5. The exact dynamic-v2.1 frozen input layer is also embedded/materialized, including the six publication-critical tables used by the dynamic pipeline.
6. Characterization, dynamic, endpoint and transition validation stages execute from the publication tree with repository-relative paths. The transition-result chain validates across v1.0, v1.1, v1.2, v1.2.1 and v1.3 frozen/checkpoint outputs.
7. The current manuscript result map uses the endpoint/dynamic chain, transition staging/factorization, KL convergence and nonequilibrium/Fisher analyses. Historical pathwise-temporal-order and PI/Fourier packages are not linked to any current manuscript result and are therefore outside the publication-critical DAG.

## Stage status

| ID | Stage | Final source/version | Status | Publication action |
|---|---|---|---|---|
| 00 | Raw current-clamp + QC/provenance | `NeuroThermo_stage1_qc_fixed` + raw ABF snapshot | COMPLETE_RELEASE_VALIDATED | clean-clone raw integrity and exact scientific QC-decision reproduction pass |
| 01 | Frozen event/QC calibration | exact server archive + canonical v3.5 names | COMPLETE_RELEASE_VALIDATED | retained with SHA-256 verification and deterministic materialization |
| 02 | Restricted HR fit | `NeuroThermo_cell_fit_v3_9_frozen_exact` | COMPLETE_INPUTS_SOURCE_AND_FULL_FROZEN_RESULTS | full expensive refit/identifiability commands are exposed; routine CI validates exact inputs and frozen outputs |
| 03 | Characterization | `NeuroThermo_characterization_v1_0` | COMPLETE_RELEASE_VALIDATED | complete v3.9 inputs bundled; compatibility XLSX generated deterministically from canonical CSV |
| 04 | Support-restricted dynamics | `NeuroThermo_dynamic_v2_1` | COMPLETE_RELEASE_VALIDATED | exact six-file frozen input layer bundled/materialized and validated |
| 05 | Endpoint uncertainty ensemble | `NeuroThermo_endpoint_ensemble_v1_0_1` | COMPLETE_RELEASE_VALIDATED | frozen layer and validation path present |
| 06 | Transition base | `NeuroThermo_transition_v1_0` | COMPLETE_RELEASE_VALIDATED | portable config/input validation present; heavy 988-scenario recompute excluded from routine CI |
| 07 | Corrected transition staging | `NeuroThermo_transition_v1_1` | COMPLETE_RELEASE_VALIDATED | zero-new-simulation reprojection and staging validation pass |
| 08 | Intrinsic × drive transition | `NeuroThermo_transition_v1_2` + `v1_2_1` | COMPLETE_RELEASE_VALIDATED | v1.1→v1.2 assembly and frozen/result validation pass; heavy surface recompute checkpoint/resume only |
| 09 | Factorized drive | `NeuroThermo_transition_v1_3_frozen_exact` | COMPLETE_RELEASE_VALIDATED | exact source/checkpoint chain validates across all 988 v1.2 scenarios |
| 10 | Full-coverage KL | `NeuroThermo_KL_convergence_v1_0_1` | COMPLETE_RELEASE | source/config/results bundled; final clean-clone full-run evidence still required on release SHA |
| 11 | Nonequilibrium geometry/Fisher | `NeuroThermo_nonequilibrium_geometry_v1_0_1` | COMPLETE_RELEASE | source/config/results bundled; final clean-clone full-run evidence still required on release SHA |
| 12 | Figure source-data assembly | `code/assemble_figure_source.py` | COMPLETE_RELEASE_VALIDATED | all renderer-consumed Fig1–3 source tables reconstruct against frozen references |
| 13 | Figure rendering | canonical Python renderers | COMPLETE_RELEASE_VALIDATED | Python/matplotlib is canonical; R retained as optional alternative renderer |
| 14 | Master workflow | `code/run_full_analyses.sh` + `code/reproduce_figures.sh` | COMPLETE_RELEASE | publication stages are wired; routine CI remains intentionally lighter than full expensive DAG |
| 15 | Environment lock | Python 3.9.25 historical environment | PYTHON_CANONICAL_COMPLETE_R_OPTIONAL | exact tested Python versions pinned; R not required for canonical reviewer path |
| 16 | Pathwise temporal order | historical auxiliary package | OUT_OF_CURRENT_MANUSCRIPT_SCOPE | do not add unless a manuscript result explicitly depends on it |
| 17 | 2D KL auxiliary | covered by KL v1.0.1 marginal analysis | COVERED_BY_CURRENT_KL_PIPELINE | no separate publication-critical package required for current result map |
| 18 | PI/Fourier | historical auxiliary package | OUT_OF_CURRENT_MANUSCRIPT_SCOPE | do not add unless manuscript scope changes |

Machine-readable details are in `PUBLICATION_WORKFLOW_INVENTORY.tsv`.

## Confirmed computational chain

```text
50 raw ABFs + QC/provenance
  -> stage1 candidate extraction / classifier / fixed visual QC
  -> exact frozen calibration and v3.5 accepted selections
  -> exact cell_fit_v3_9 source + complete frozen v3.9 outputs
  -> characterization_v1_0
  -> exact dynamic_v2_1 frozen input layer
  -> dynamic_v2_1
  -> endpoint_ensemble_v1_0_1
  -> transition_v1_0
  -> transition_v1_1
  -> transition_v1_2
  -> transition_v1_2_1
  -> transition_v1_3
  -> Fig1–3/S1 source-data assembly and rendering

final transition/support-state scenarios
  -> KL_convergence_v1_0_1
  -> current Fig4 + KL supplementary outputs

final transition/support-state scenarios
  -> nonequilibrium_geometry_v1_0_1
  -> Fisher/time-reversal/stationarity supporting outputs
```

The constructed path coordinate `p` is a coordinate through model-state space and is not interpreted as disease time, causal progression or evidence of irreversible one-way evolution.

## Key frozen/calibration hashes

- calibration archive: `0c930506021826aec8ee2987fe83cd4a1537fa42b6d3fad335a5520fcbb610bd`
- candidate events: `af35c327b313482f534aa59669a47e52a4078f912a5e342efcfddf0158455640`
- accepted sweeps v3.5: `dad46b831eb4613af4a49673f83854e4ef48b81d0934c087234562d81a447a54`
- peak overrides v3.5: `64e35808199e6108355b015b4ca9ded6070deed852927877e705ccf118e95069`
- threshold brackets v3.5: `47ba271e6b8d70704de1c49aaac3677c6ee21e3001f33faaafad8761177f9741`
- v3.1 cell summary: `85b1fa2c457e4affc0db438cf885b4406f61b943cbc08073fbdeb7f4b57f42f9`
- v3.1 sweep summary: `5663d59c35aeb105ee45b0c4c8606375210294f377a6ee3adcd771356a70ab12`
- v3.1 identifiability: `16e810e3331a0f6eb6bc1c815bb0e0d5574ee93966b95c00346522d5470957d1`
- v3.9 seed summary: `cb74bc0783c9fd1db11cacba13ccabd273cfc225e6dc019ab6e4215433dceb72`

## Remaining release blockers

The scientific/raw-to-QC portability blocker is resolved. The remaining blockers are release-engineering issues:

1. **CI coverage of release-critical changes.** The publication preflight workflow must run on every change that can alter the release snapshot, not only the current subset of paths.
2. **Whole-release manifest freshness.** `MANIFEST.sha256` must be regenerated after all code/document/workflow changes and checked on the same final SHA.
3. **Final clean-clone evidence.** The final release SHA needs a successful publication preflight; a separate full downstream clean-clone run should execute the complete non-optimizer DAG through KL, nonequilibrium and figure rendering. Heavy HR optimization/surface recomputation can remain separate because exact frozen checkpoints are part of the release.
4. **Immutable/public snapshot.** The validated release must be exposed as an immutable tag/release or separate public repository snapshot. The current source repository is private, so a publication branch inside it is not independently public.
5. **Citation metadata.** Add top-level release citation metadata (for example `CITATION.cff`) when the public repository identity/DOI is finalized.

## Next remediation order

1. Broaden publication CI triggers and add release-manifest checking to the final gate.
2. Run the complete clean-clone downstream replay/validation chain on the resulting exact SHA.
3. Regenerate `MANIFEST.sha256` last and re-run all final gates on that manifest-bearing SHA.
4. Create the immutable public snapshot/tag/release and citation metadata.
