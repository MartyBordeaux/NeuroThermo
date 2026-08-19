# NeuroThermo dependency-aware fingerprint v0.5.1

## Analysis unit

- Recorded cells are independent observations: 13 WT and 7 SCA3.
- Animal identity is removed at the input boundary and is not used.
- Scope of inference: the recorded cells, not the animal population.

## Dependency-reduced primary fingerprint

- Median burden difference (SCA3 - WT): 2.46746.
- Exact recomputed two-sided permutation p: 0.000193498 across 77520 labelings.
- In-sample WT-referenced AUC: 1.0000.
- Default WT-exit calls: 0 WT and 4 SCA3 cells.

## Dependency audit

- Absolute rheobase median difference: 0 pA; exact p=1.
- Capacitance-normalized rheobase difference: 1.77789 pA/pF; exact p=0.00469556. It is a derived diagnostic, not independent core evidence.
- Work per spike difference: 1503.4 fJ; exact p=0.0731682. It remains a derived diagnostic linked to spike count.

## Robustness

- Specification AUC range: 0.9670 to 1.0000.
- Current-window exact p range: 0.000116099 to 0.000348297.
- Leave-one-cell AUC range: 1.0000 to 1.0000.

The endpoint q coordinate and WT-exit boundary are operational phenotype summaries. They are not disease probabilities, biological time, or evidence of irreversible degeneration.
