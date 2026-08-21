# NeuroThermo v0.7.1 — surrogate validation of ordinal predictive information

## Frozen design

The analysis unit is one cell. The frozen cohort contains all 13 WT and 7 SCA3 cells. No animal-level inference is performed. The current grid is 100–600 pA. The primary null is value shuffling, which quantifies finite-sample and overlapping-window bias of ordinal predictive information. AAFT is a secondary null preserving the amplitude distribution and approximately the power spectrum.

## Primary endpoint

Primary metric: observed ordinal PI minus the median shuffled-surrogate PI. WT-minus-SCA3 cell-AUC difference: 0.046173299; exact p=1.2899897e-05.

## Secondary checks

AAFT-excess WT-minus-SCA3 cell-AUC difference: 0.034015358; exact p=9.0299278e-05.

Surrogate-corrected adjusted burden curve AUC difference: 1.2581197; exact p=0.00042569659.

Secondary adjusted exits: SCA3_02=500 pA.

## Reproduction requirement

PI is recomputed from the raw stationary ABF interval with the frozen v0.3.1 estimator. The maximum permitted absolute discrepancy from v0.3.1 is 1e-09 nats. Failure of this check aborts the run.

## Interpretation boundary

The primary result tests whether the WT–SCA3 difference survives correction for estimator bias caused by finite record length and overlapping ordinal windows. AAFT sensitivity asks whether the difference remains beyond structure captured by the amplitude distribution and approximately by the spectrum. Neither result establishes disease time, irreversibility, causal mechanism, or a thermodynamic phase transition. I_exit remains a secondary current-stress threshold and is not a disease-onset time.
