# Transition ensemble v1.2.1

No new Hindmarsh–Rose simulations are performed. The frozen v1.2 checkpoint surfaces are reprocessed with scenario-first uncertainty propagation.

Correct order: scenario crossing -> v1.1-compatible within-pair weighted marker distribution with crossing-support weight retained separately -> equal-weight biological-pair ensemble.

## Recovery of frozen v1.1 one-dimensional paths

- coupled / WT_exit: surface=0.391909, frozen v1.1=0.398286, |delta|=0.006377.
- coupled / balance: surface=0.662635, frozen v1.1=0.676380, |delta|=0.013746.
- coupled / SCA3_entry: surface=0.837546, frozen v1.1=0.837017, |delta|=0.000529.
- drive_early / WT_exit: surface=0.233456, frozen v1.1=0.233544, |delta|=0.000087.
- drive_early / balance: surface=0.522213, frozen v1.1=0.522569, |delta|=0.000356.
- drive_early / SCA3_entry: surface=0.716332, frozen v1.1=0.716867, |delta|=0.000535.
- drive_late / WT_exit: surface=0.289262, frozen v1.1=0.400945, |delta|=0.111684.
- drive_late / balance: surface=0.765491, frozen v1.1=0.793264, |delta|=0.027773.
- drive_late / SCA3_entry: surface=0.873497, frozen v1.1=0.878271, |delta|=0.004775.

## Scenario-aware drive sensitivity at ISI stage boundaries

- SCA3_entry: drive dominance median=0.5180; |dA/ddrive|=1.3034; |dA/dintrinsic|=0.9893.
- WT_exit: drive dominance median=0.3389; |dA/ddrive|=0.7692; |dA/dintrinsic|=0.8161.
- balance: drive dominance median=0.4028; |dA/ddrive|=0.9994; |dA/dintrinsic|=0.9342.

Parameter-support states without a marker are represented by the separately reported crossing-support weight; marker quantiles use the same finite-marker weighting convention as frozen v1.1 for direct comparability.
