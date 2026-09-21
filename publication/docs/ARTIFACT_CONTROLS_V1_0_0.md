# NeuroThermo artifact controls v1.0.0

This is the prespecified reviewer-facing falsification audit added after the final route-robustness and diffusion-sensitivity analyses.

## Scope

The pipeline tests whether the state-before-firing ordering survives:
- pair-specific endpoint re-anchoring;
- Jensen-Shannon, Hellinger, total-variation, kNN-KL, and sliced-Wasserstein-2 state distances;
- all 28 WT-WT and 6 SCA3-SCA3 unordered null pairs, both orientations;
- sample-window, histogram-bin, pseudocount, and state-marginal sensitivity;
- slow-variable stationarity and a 12-s long-burn audit;
- leave-one-cell-out analysis over the 12 core cells;
- unsupported axis-crossing classification;
- endpoint residual diagnostics;
- secondary fixed-J sensitivity.

The frozen route ensemble is the same 525-route family used by the final route-robustness analysis.

## Frozen decision rule

The only global verdicts allowed before calculation were:
- `ARTIFACT_CONTROLS_SUPPORT_PRIMARY_ORDERING`;
- `ARTIFACT_SENSITIVITY_OR_UNRESOLVED`.

No failed gate was eligible for post-hoc threshold or estimator retuning.

## Result

Final verdict: `ARTIFACT_SENSITIVITY_OR_UNRESOLVED`.

Passed components:
- endpoint re-anchoring;
- symmetric metrics;
- state marginals;
- sampling convergence;
- histogram sensitivity;
- non-histogram cross-checks;
- leave-one-cell-out;
- long-burn stationarity;
- half-window stationarity.

Unresolved component:
- same-genotype nulls. The SCA3-SCA3 null remained near zero, whereas WT-WT showed non-zero within-WT heterogeneity and failed the prespecified near-zero gate. The WT-WT shift did not reproduce the negative WT-to-SCA3 state-before-firing ordering.

## Reproducibility

Source archive:
`publication/code_archives/neurothermo_artifact_controls_v1_0_0_source_only.zip`

SHA-256:
`006d346ed060da9cb01f96448f8e5c6494ea53f2e441bf05e096a0f170e122a7`

Complete calculation package SHA-256:
`8715ab05bd6efd701d2e4655b2fd444b0d500e1a39e6a1547188edeed7c9a9fe`

Frozen 525-route table SHA-256:
`132d85c9fe10c5b9e50e2cbf0e9ebcb816c084deadf88e8841ec31f26e470d5e`

The source-only archive omits the large frozen route binary and smoke binary input only to avoid duplicating large immutable artifacts; it retains their expected paths/checksums and all executable source, configurations, tests, documentation, and launch scripts.
