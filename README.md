# NeuroThermo

Public reproducibility repository for the manuscript:

**Dynamical state reorganization precedes the firing phenotype in WT-SCA3 Purkinje-cell model space**

## Publication package

The current manuscript-facing entry point is:

**[`publication/`](publication/)**

It contains the final result-to-code map, checksums, source-code archive for the four manuscript-facing post-rejection pipelines, and compact diffusion-sensitivity results.

The repository also retains the frozen upstream analysis chain: raw current-clamp recordings, QC/provenance material, frozen spike selections, restricted four-parameter Hindmarsh-Rose fits, endpoint ensembles, transition analyses, and supporting result trees.

## Current scientific result

Within the endpoint-constrained intrinsic-drive Hindmarsh-Rose model space, the full-state KL balance boundary is encountered before the reduced firing-phenotype balance boundary across a broad prespecified family of admissible monotone routes. The final robustness analysis uses 525 routes, including 512 random monotone routes, and the ordering remains stable when all three stochastic diffusion coefficients are halved or doubled.

Model-space coordinates are not disease time. Cross-combinations, scenarios, stochastic seeds, atlas nodes, and routes are computational units rather than independent biological replicates.

## Frozen model

- fitted parameters: `b`, `r`, `s`, `kappa_I`;
- fixed constants: `a=1`, `c=1`, `d=5`, `x_R=-1.6`;
- primary multi-sweep cohort: 18 cells (12 WT, 6 SCA3);
- core model-space subset: 8 WT + 4 SCA3 cells;
- 32 dependent WT x SCA3 endpoint combinations;
- 264 retained support scenarios;
- 121-node intrinsic x drive atlas;
- 525 frozen routes.

Historical directory and branch names are retained to preserve exact provenance of earlier frozen analyses. New manuscript-facing material is organized only under the neutral `publication/` entry point.
