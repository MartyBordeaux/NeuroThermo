# Endpoint ensemble v1.0.1

## Purpose

Endpoint ensemble v1.0.1 converts frozen v3.9 + dynamic v2.1 results into uncertainty-aware WT and SCA3 endpoint clouds for the next transition/staging analysis.

No additional fitting is performed.

## Biological and model support

- Primary biological units: 18 multi-sweep cells (12 WT, 6 SCA3).
- Actual HR support states: 64 total.
- Cells with at least one near-optimal HR alternative: 13/18.
- Near-optimal solutions are a within-cell model-uncertainty layer, not additional biological replicates and not a Bayesian posterior.
- Each cell has equal biological weight within its genotype group; that cell's weight is divided among its best and near-optimal support states.

## Full observable record

Each cell/support state stores:

- capacitance-normalized rheobase `J_rheo`;
- firing rate at experiment-supported `q=0.75`;
- mean ISI at `q=0.75`;
- firing rate at `q=0.50` when supported;
- mean ISI at `q=0.50` when supported;
- HR coordinates `(b,r,s,kappa_I)`.

No `q=0.50` imputation is used.

## Primary transition projection

The primary two-dimensional transition space is

`(log10 J_rheo, log10 firing_rate_q75)`.

Coordinates are robust-scaled using the pooled best-cell median and MAD. Mean ISI at `q=0.75` is retained as a parallel projection/validation observable because transient spike trains can make firing rate and local ISI non-equivalent.

The experimental best-cell WT/SCA3 centroid distance is

`D = 1.831 pooled-MAD units`.

Best-fit model centroid distance is `1.798`, and the best-model transition direction differs from the experimental direction by only about `1.4 degrees`.

## Endpoint separation

Experimental best-cell medians:

- `J_rheo`: SCA3 `2.650`, WT `0.642 pA/pF`;
- `f_q75`: SCA3 `34.22`, WT `75.33 Hz`;
- `ISI_q75`: SCA3 `25.93`, WT `12.62 ms`.

The deterministic sensitivity envelope over the available within-cell near-optimal support states does not reverse the primary group direction:

- SCA3-WT standardized `log10 J_rheo` difference remains in `[+1.386, +1.636]`;
- SCA3-WT standardized `log10 f_q75` difference remains in `[-2.447, -0.644]`.

These are sensitivity envelopes over enumerated support states, not confidence or credible intervals.

## Uncertainty relative to biological spread

Median within-cell support span / between-cell best-solution IQR on the log scale:

- rheobase: `0.134` WT, `0.129` SCA3;
- `q=0.75` firing rate: `0.363` WT, `0.348` SCA3;
- `q=0.75` mean ISI: `0.374` WT, `0.444` SCA3;
- `q=0.50` firing rate: `0.217` WT, `3.326` SCA3.

This is why `q=0.75` is the primary transition anchor and `q=0.50` is retained only as an extended support layer.

## Secure subset

At the adopted 20% diagnostic robustness threshold:

- rheobase secure: 15/18;
- `q=0.75` firing rate secure: 14/18;
- `q=0.75` mean ISI secure: 14/18;
- joint core secure: 12/18 (8 WT, 4 SCA3).

No cell is excluded from the full endpoint ensemble because of model uncertainty; uncertain cells simply carry a broader set of admissible support states.

## Transition readiness

The next WT→SCA3 pipeline should use biological pair weighting first and divide that pair weight across parameter-support scenarios. It should preserve two parallel observable projections:

1. `(log10 J_rheo, log10 firing_rate_q75)`;
2. `(log10 J_rheo, log10 mean_ISI_q75)`.

Staging conclusions should be considered stronger when reproduced in both projections.
