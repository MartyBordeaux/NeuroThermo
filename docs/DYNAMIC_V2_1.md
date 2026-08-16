# Dynamic characterization v2.1

## Purpose

Dynamic v2.1 is a post-fit analysis on frozen v3.9 solutions. It was introduced after v2.0 showed that fixed-`lambda = J/J_rheo` comparisons could fall outside the experimentally sampled current range and could mix rheobase uncertainty with suprathreshold-dynamics uncertainty.

No HR parameters are refit in v2.1.

## Experimental-support coordinate

For each primary multi-sweep cell,

`q = (J - J_rheo,best) / (J_max,obs - J_rheo,best)`.

Only interpolation between experimentally observed spiking currents is allowed. Extrapolation is prohibited.

Support coverage:

- `q=0.25`: 12/18 cells (3/6 SCA3, 9/12 WT);
- `q=0.50`: 17/18 cells (5/6 SCA3, 12/12 WT);
- `q=0.75`: 18/18 cells (6/6 SCA3, 12/12 WT).

Near-optimal solutions are evaluated at the same physical current as the corresponding best-fit solution; an alternative solution's own rheobase does not silently change the applied current in the suprathreshold robustness test.

## Primary experiment-supported results

At `q=0.50`, recorded-cell medians are:

- firing rate: SCA3 `28.85 Hz`, WT `55.24 Hz`;
- mean ISI: SCA3 `31.02 ms`, WT `15.94 ms`.

At `q=0.75`, recorded-cell medians are:

- firing rate: SCA3 `34.22 Hz`, WT `75.33 Hz`;
- mean ISI: SCA3 `25.93 ms`, WT `12.62 ms`.

Best-fit model medians at `q=0.75` preserve the same direction:

- firing rate: SCA3 `42.58 Hz`, WT `75.20 Hz`;
- mean ISI: SCA3 `22.05 ms`, WT `12.28 ms`.

The best-fit rheobase medians are:

- SCA3 `2.650 pA/pF`;
- WT `0.642 pA/pF`.

Thus the supported phenotype separates into two axes: reduced absolute excitability (higher capacitance-normalized rheobase) and slower suprathreshold spiking dynamics (lower firing rate / longer ISI).

## Robustness to near-optimal HR solutions

At `q=0.75`, 13 cells have evaluable near-optimal alternatives. Under the diagnostic 20% robustness threshold:

- firing rate is stable in 9/13 evaluable cells;
- mean ISI is stable in 9/13;
- train duration is stable in 11/13.

Five cells have no supported near-optimal alternative at this point because their v3.9 solution is practically identifiable under the adopted alternative-search criterion.

When practical parameter identifiability and near-optimal robustness are combined, endpoint ensemble v1.0.1 classifies 14/18 cells as secure for `q=0.75` firing rate, 14/18 for mean ISI, and 12/18 for the joint rheobase+rate+ISI core.

## Phase-resolved status

Phase-resolved HR geometry is retained as exploratory rather than primary. Among evaluable near-optimal solutions across supported q-points, the 20% stability rate is substantially lower for full phase profiles than for the primary scalar phenotype; in particular, absolute `z(phi)` is not robust. Therefore phase-profile claims should not be used as the principal biological endpoint before uncertainty is propagated explicitly.

## Statistical scope

Animal identifiers remain incomplete for part of WT and only two recovered animals per group are available. Formal animal-level genotype p-values are not reported. Cell-level effects describe the recorded-cell ensemble; recovered-animal summaries are descriptive.
