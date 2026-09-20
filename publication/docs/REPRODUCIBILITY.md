# Reproducibility notes

## Frozen analysis source

The four final manuscript-facing post-rejection pipelines are preserved in one verified source archive:

`publication/code_archives/publication_source_code.tar.xz`

Because the update interface accepted text files only, the binary archive is stored as seven ordered Base64 parts. Running `publication/code_archives/reconstruct_source_code.sh` concatenates the parts, decodes the archive, verifies SHA-256, and extracts the source tree.

Archive SHA-256:

`8a2b142691af14f528d13847371e8eeed5c002577d1253fa67022f6fa4b46241`

The archive contains the exact source/configuration/documentation material recovered from:
- `neurothermo_route_robustness_v1_0_0`
- `neurothermo_current_intervention_predictions_v1_0_3`
- `neurothermo_raw_proxy_robustness_v1_0_0`
- `neurothermo_noise_landscape_robustness_v1_0_0`

## Source archive fingerprints

- route robustness: `3db44bca3f7cb2a56f099390d6acdeeb41d4e65c2ae82595593824aaf5a0cb8e`
- current intervention: `e80cac8153b2e44e2d9a4fff0f5163c0fac98fcf7fa9eee5113ef346b764ce58`
- raw proxy robustness: `d7cb3689738be4d2c07f6fb01c4cdbd4e296da6dd6555836fe692d7c024a349a`
- noise landscape robustness: `84d89b2b732a69e3fdc85c16658a524f9322aedd512e3c5b9b86d9050f2c8008`

## Noise-landscape result fingerprint

Compact result archive SHA-256:

`b4419665af206a71db4bd30ddb58fdda26bcc5806ace4d9e2545a317f1b70ef3`

The key manuscript-facing tables from this archive are committed as plain text under `publication/results/noise_landscape_robustness_v1_0_0/`.

## Large result archives

The full route-robustness/current-intervention/raw-proxy result trees are larger binary artifacts retained in the project archive. The public Git repository contains their source code, frozen upstream data/results, and manuscript-facing compact outputs/provenance. A future immutable archival release can additionally mirror the complete binary result trees without changing the analysis identity.
