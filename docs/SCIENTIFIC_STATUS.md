# Scientific status after v3.9, dynamic v2.1, endpoint ensemble v1.0.1, and transition v1.3.1

Primary multi-sweep cohort: 18 cells, comprising 12 WT and 6 SCA3. Two additional accepted single-sweep cells are retained only as secondary descriptive records.

## Final cell-fit status

The stable group-associated HR coordinates remain increased `s` and a more moderate increase in `b` in the recorded SCA3 cells. The apparent WT–SCA3 difference in `r` does not survive wide-bound stress and `r` remains the least identifiable coordinate. Raw `kappa_I` strongly tracks capacitance and is treated as a model input-scaling coordinate rather than an independently established SCA3 biophysical phenotype.

## Primary robust dynamical phenotype

Dynamic v2.1 restricts analysis to experimentally supported suprathreshold currents and does not extrapolate beyond the observed spiking-current range.

At `q=0.75`, supported in all 18 primary cells:

- firing rate: SCA3 median `34.22 Hz`, WT `75.33 Hz`;
- mean ISI: SCA3 `25.93 ms`, WT `12.62 ms`.

Capacitance-normalized best-fit rheobase medians are `2.650 pA/pF` in SCA3 and `0.642 pA/pF` in WT.

The working phenotype is therefore two-dimensional: reduced absolute excitability together with slower suprathreshold spike-train dynamics.

## Robustness to parameter non-identifiability

Endpoint ensemble v1.0.1 contains 64 admissible HR support states from 18 biological cells. Near-optimal parameter solutions are a within-cell uncertainty layer and do not receive independent biological weight. The primary endpoint direction remains stable across the enumerated support states.

## Frozen WT→SCA3 staging geometry

Transition v1.1 corrected the secondary rate projection and retained latency-invariant ISI as the primary staging observable. In the core-secure reference set:

- WT-exit: `A_ISI = 0.135829`;
- balance: `A_ISI = 0.500000`;
- SCA3-entry: `A_ISI = 0.797856`.

Across the 32 core-secure WT×SCA3 biological pairs, primary ISI staging medians are:

| path family | WT-exit | balance | SCA3-entry |
|---|---:|---:|---:|
| drive early | 0.234 | 0.523 | 0.717 |
| coupled | 0.398 | 0.676 | 0.837 |
| drive late | 0.401 | 0.793 | 0.878 |

Thus the ordering WT-exit < balance < SCA3-entry is preserved, while the timing of the input-related transition shifts the macroscopic transition position.

## Scenario-first intrinsic × drive map

Transition v1.2.1 corrected uncertainty aggregation by finding crossings at the individual support-state scenario level before within-pair and across-pair aggregation.

Near the three ISI boundaries in the core-secure subset, median local sensitivities are:

| stage | |dA/d drive| | |dA/d intrinsic| | drive-dominance fraction |
|---|---:|---:|---:|
| WT-exit | 0.769 | 0.816 | 0.339 |
| balance | 0.999 | 0.934 | 0.403 |
| SCA3-entry | 1.303 | 0.989 | 0.518 |

The transition is therefore not uniformly drive-dominated. Intrinsic and drive effects are comparable early and around balance, whereas relative drive sensitivity increases toward SCA3-entry.

The 31×31 2D grid reproduces the frozen v1.1 early/coupled paths closely except for `drive_late WT-exit`, where the absolute discrepancy is `0.112`. This is retained as a grid-resolution diagnostic; v1.1 remains authoritative for that particular one-dimensional early crossing.

## Factorized drive decomposition

Transition v1.3 decomposes the combined drive into fitted `kappa_I` and experimenter-imposed current `J`.

For core-secure ISI staging bands, the median absolute component effects are:

| stage | |ΔA_kappa| | |ΔA_J| | |interaction| | signed interaction |
|---|---:|---:|---:|---:|
| WT-exit | 0.309 | 0.215 | 0.032 | -0.040 |
| balance | 0.295 | 0.152 | 0.031 | -0.025 |
| SCA3-entry | 0.278 | 0.151 | 0.033 | -0.022 |

The `kappa_I` contribution is larger than the `J` contribution at all three primary stage boundaries. The interaction is modest and negative in these boundary regions, indicating weak sub-additivity rather than strong positive synergy.

Along the coupled trajectory, increasing `kappa_I` progress acts in the SCA3-directed sense whereas changing `J` toward the SCA3 endpoint generally acts in the opposite direction. Consequently, the experimental current protocol partially compensates the fitted input-scaling change rather than generating the staged phenotype by itself.

`J` must not be interpreted as a disease parameter: it is imposed experimentally. `kappa_I` must also not be equated directly with a membrane conductance or channel property; it is a fitted HR input-scaling coordinate and is strongly associated with capacitance.

## Figure correction v1.3.1

The numerical v1.3 results are unchanged. The three boundary plots in v1.3.0 were blank because the plotting code used pandas attribute access `b.mode`, which resolves to the `DataFrame.mode()` method rather than the `mode` column. v1.3.1 fixes only this visualization bug by using explicit column indexing and regenerates the boundary figures.

## Phase-resolved status

Best-fit HR trajectories show phase-localized differences, but full phase profiles remain substantially less robust to near-optimal HR alternatives than rheobase and spike-train timing. Absolute `z(phi)` is particularly unstable. Phase-resolved geometry remains supporting/exploratory unless uncertainty is propagated explicitly.

## Statistical scope

Only two recovered animals are currently available per group, and several WT animal IDs remain unresolved. Formal animal-level genotype p-values are not reported. Cell-level effects describe the recorded-cell ensemble; recovered-animal summaries remain descriptive.

## Next analysis boundary

The uncertainty-aware transition geometry and drive decomposition are now frozen. The next major analysis can overlay thermodynamic quantities (Fisher information, entropy/KL structure, EPR or related irreversibility measures) on the coupled transition while using the intrinsic–`kappa_I` map as the main mechanistic sensitivity layer.
