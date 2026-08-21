# Raw-data thermodynamic phenotype track through v0.7.3

## Scope

This track analyzes the recorded current-clamp traces directly, without fitting a Hindmarsh–Rose model. It is complementary to the model-based HR ensemble analysis already frozen in this repository.

The frozen cohort contains all 20 accepted recordings: 13 WT and 7 SCA3 cells. The recorded cell is the independent analysis unit. Animal-level inference is absent. The imposed current `I_ext` is a stress coordinate controlled by the experimenter and is not treated as a disease phenotype.

## Frozen milestones

### v0.3.1 — current-resolved phenotype extraction

The pipeline converts curated ABF sweeps into cell- and current-resolved scalar, integrated and dynamic phenotypes on the common 100–600 pA grid. It freezes provenance, cell identity, current support, sweep-level quality control and the base phenotype tables used by later stages.

### v0.5.1 — dependency-reduced endpoint fingerprint

The endpoint fingerprint uses one representative feature per physiological domain to reduce deterministic and near-deterministic dependencies. The median composite burden difference was 2.46746 robust z units in the SCA3 direction. The exact two-sided cell-label permutation p-value was 0.000193498 across all 77,520 labelings. In-sample WT-referenced AUC was 1.0; this is an endpoint separation statistic, not external classifier validation.

### v0.6.1 — current-resolved vulnerability

Firing output and conditional timing were combined into one excitability/timing domain, while predictive dynamics remained a second domain. A strict stress exit required both domains and their mean burden to exceed the frozen threshold, a target-specific leave-one-out WT rank criterion, and persistence across adjacent current levels.

The burden-curve AUC difference was 1.40277 in the SCA3 direction with exact p=1.28999e-5. No WT cell and two SCA3 cells crossed the calibrated boundary by 600 pA. The restricted-mean `I_exit` difference was 50 pA with exact p=0.00130289. `I_exit` is a threshold under imposed current stress, not degeneration time.

### v0.7.0 — activity- and technical-adjusted predictive dynamics

Predictive information was residualized at each current by label-blind leave-one-cell-out ridge regression using firing rate, mean ISI and its missingness indicator, baseline noise and stationary sample count. The adjusted two-domain burden AUC difference remained significant (exact p=0.000386997), and the residual predictive-dynamics AUC difference had exact p=0.00216718. This supports a temporal-organization contribution not explained by the included activity and technical covariates.

### v0.7.1 — shuffled ordinal predictive-information null

The adjacent-window ordinal PI analysis produced a raw WT-minus-SCA3 cell-AUC difference of 0.0462058 nats and a shuffle-centered difference of 0.0461733 nats. Both had exact expected-direction p=1.28999e-5. Because adjacent ordinal words shared three of four samples, the analysis was treated as preliminary and replaced by non-overlapping words in v0.7.2.

### v0.7.2 — non-overlapping PI and failed IAAFT sensitivity

The primary non-overlapping code-pair lag was 4 ms; 8, 16 and 32 ms were secondary. The frozen shuffle-centered WT-minus-SCA3 cell-AUC differences were:

| lag | difference (nats) | exact directional p | maxT p across lags |
|---:|---:|---:|---:|
| 4 ms | 0.14784769 | 0.00918473 | 0.00936533 |
| 8 ms | 0.13996580 | 0.00666925 | 0.01403509 |
| 16 ms | 0.04761938 | 0.03397833 | 0.26495098 |
| 32 ms | 0.02864381 | 0.02497420 | 0.36842105 |

The IAAFT amplitude and autocorrelation checks passed, but the spectral-amplitude NRMSE and convergence gates failed. IAAFT results are excluded from inference. The failure is preserved as provenance because it motivated the exact-spectrum construction in v0.7.3.

### v0.7.3 — rank-Gaussian exact-spectrum Fourier validation

Each stationary voltage trace was transformed monotonically to empirical normal scores. This leaves ordinal codes invariant. Random Fourier phases then generated surrogates preserving the complete rFFT magnitude spectrum of the rank-Gaussian trace.

Ordinal invariance held to 2.22e-16 nats. Maximum spectral-amplitude NRMSE was 4.14e-16 and maximum circular-ACF error was 5.55e-16. All exact-spectrum fidelity gates passed.

| lag | Fourier-centered WT−SCA3 AUC difference (nats) | exact directional p | maxT p across lags |
|---:|---:|---:|---:|
| 4 ms | 0.13650280 | 0.00959752 | 0.00959752 |
| 8 ms | 0.12458749 | 0.00764964 | 0.01733746 |
| 16 ms | 0.04288244 | 0.04369195 | 0.26724716 |
| 32 ms | 0.02766009 | 0.02765738 | 0.35912023 |

After maxT correction across four lags and eleven currents, the exact-spectrum effect remained at 450 and 500 pA for the 4-ms phenotype and at 450 pA for the 8-ms phenotype. Leave-one-cell-out 4- and 8-ms effects were positive for every omitted cell. The 4-ms cell-AUC pairwise ordering probability was 75/91=0.824, so the group phenotype is strong but not a perfect cell classifier.

## Scientific interpretation

WT cells show more short-timescale ordinal temporal structure than SCA3 cells after subtracting structure reproduced by both a shuffle null and an exact-spectrum linear Gaussian-copula null. The group difference is concentrated at 4–8 ms and is revealed most strongly under 450–500 pA current stress.

This result identifies higher-order temporal organization not reproduced by the frozen null models. It does not establish entropy production, thermodynamic irreversibility, a causal degeneration mechanism, physical disease time or a phase transition.

## Frozen inference boundary

- Unit of inference: recorded cell.
- Cohort: 13 WT and 7 SCA3 cells; no cell is discarded for low group size.
- Exact labelings: 77,520.
- Current is an imposed perturbation coordinate.
- IAAFT is excluded from inference after failed fidelity gates.
- v0.7.3 Fourier-centered non-overlapping PI at 4 ms is the frozen primary temporal-order phenotype.

