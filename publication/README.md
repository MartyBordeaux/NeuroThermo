# NeuroThermo publication package

This directory is the manuscript-facing entry point for the study:

**Dynamical state reorganization precedes the firing phenotype in WT-SCA3 Purkinje-cell model space**

The repository is public and contains the experimental current-clamp data, frozen QC and spike-selection inputs, Hindmarsh-Rose fitting outputs, endpoint/transition results, final post-rejection analysis source code, and manuscript-facing result tables.

## Final analyses used by the manuscript

- `neurothermo_route_robustness_v1_0_0` - 121-node intrinsic x drive landscape, 525-route robustness analysis, contour ordering, discordance geometry, and supporting factorial attribution.
- `neurothermo_current_intervention_predictions_v1_0_3` - physical-current predictions with SCA3 endpoint parameters fixed.
- `neurothermo_raw_proxy_robustness_v1_0_0` - retrospective raw-voltage proxy robustness audit.
- `neurothermo_noise_landscape_robustness_v1_0_0` - confirmatory diffusion-coefficient sensitivity on the final 2D landscape.
- `neurothermo_artifact_controls_v1_0_0` - prespecified reviewer-facing artifact controls: endpoint re-anchoring, alternative state distances, same-genotype nulls, density/sampling sensitivity, stationarity, and leave-one-cell-out.

The source files for these four final pipelines are frozen in `code_archives/publication_source_code.tar.xz`, stored losslessly as Base64 parts because the repository connector used for this update accepts text files only. Reconstruct and verify the archive with:

```bash
cd publication/code_archives
bash reconstruct_source_code.sh
```

The reconstructed archive contains 137 source/configuration/documentation files and has SHA-256:

`8a2b142691af14f528d13847371e8eeed5c002577d1253fa67022f6fa4b46241`

## Manuscript-facing results

The final diffusion-sensitivity tables are under:

`results/noise_landscape_robustness_v1_0_0/`

The primary verdict is:

`DIFFUSION_INTENSITY_ROBUST_ON_2D_LANDSCAPE`

Halving or doubling all three diffusion coefficients preserved the state-before-firing sign for all 525 tested routes; KL-field sign agreement with baseline exceeded 98%, and the reverse ensemble discordance region remained zero.

The earlier route-robustness, current-intervention, and raw-proxy result trees are already present in the repository's frozen analysis history. The result-to-code mapping in `docs/RESULT_TO_CODE_MAP.md` identifies the authoritative pipeline for each manuscript claim.

## Upstream frozen material

The repository contains the original current-clamp ABF data, QC/provenance material, frozen accepted spiking sweeps, peak overrides, threshold brackets, the restricted four-parameter Hindmarsh-Rose fits, and frozen endpoint/transition analyses. Historical directory names are retained to preserve reproducibility, but this `publication/` directory is the current neutral manuscript-facing entry point.

## Interpretation boundary

The principal conclusion is conditional on the endpoint-constrained intrinsic-drive model space and the prespecified family of admissible routes. Model-space coordinates are not disease time. Cross-combinations, scenarios, stochastic seeds, atlas nodes, and routes are computational units rather than biological replicates.

## Checksums and provenance

See:
- `CHECKSUMS.sha256`
- `docs/RESULT_TO_CODE_MAP.md`
- `docs/REPRODUCIBILITY.md`

## Artifact-control audit

The final reviewer-facing artifact-control source is archived at `code_archives/neurothermo_artifact_controls_v1_0_0_source_only.zip`. Scientific design and frozen gates are described in `docs/ARTIFACT_CONTROLS_V1_0_0.md`. The original complete package used for the calculation had SHA-256 `8715ab05bd6efd701d2e4655b2fd444b0d500e1a39e6a1547188edeed7c9a9fe`; the committed source-only archive has SHA-256 `006d346ed060da9cb01f96448f8e5c6494ea53f2e441bf05e096a0f170e122a7`.
