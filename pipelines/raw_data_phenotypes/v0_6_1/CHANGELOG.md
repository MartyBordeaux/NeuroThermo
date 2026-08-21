# Changelog

## 0.6.1

- Replaced the self-containing WT maximum with target-specific nested leave-one-out WT rank calibration.
- Combined firing rate and mean ISI into one excitability/timing domain.
- Required independent confirmation from the predictive-dynamics domain for strict exit.
- Recomputed target-specific scoring and calibration under every exact whole-cell label permutation.
- Recalculated the full WT rank calibration inside every bootstrap replicate.
- Aligned plotted and tested group curves to the same target-specific score.
- Added expected-direction feature-level p-values and output hashes.

## 0.6.0

- Added cell-level current-resolved dynamic vulnerability trajectories from 100 to 600 pA.
- Added a frozen, persistent strict `I_exit` and envelope-only diagnostic exit.
- Added exact whole-cell-trajectory permutation inference over all 77,520 labelings.
- Added current-wise maxT correction, curve-AUC tests and raw-feature diagnostic curves.
- Added WT-reference bootstrap stability for individual cells.
- Explicitly removed animal-level inference and dropped `animal_id` at input.
- Retained Python 3.9 compatibility and package-local `results/` output.
