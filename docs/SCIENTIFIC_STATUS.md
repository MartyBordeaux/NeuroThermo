# Scientific status after v3.9, dynamic v2.1, and endpoint ensemble v1.0.1

Primary multi-sweep cohort: 18 cells, comprising 12 WT and 6 SCA3. Two additional accepted single-sweep cells are retained only as secondary descriptive records.

## Final cell-fit status

The stable group-associated HR coordinates remain increased `s` and a more moderate increase in `b` in the recorded SCA3 cells. The apparent WT–SCA3 difference in `r` does not survive wide-bound stress and `r` remains the least identifiable coordinate. Raw `kappa_I` strongly tracks capacitance and is treated primarily as an input-scaling/nuisance coordinate rather than an independent SCA3 phenotype.

## Primary robust dynamical phenotype

Dynamic v2.1 restricts analysis to experimentally supported suprathreshold currents and does not extrapolate beyond the observed spiking-current range.

At `q=0.75`, supported in all 18 primary cells:

- firing rate: SCA3 median `34.22 Hz`, WT `75.33 Hz`;
- mean ISI: SCA3 `25.93 ms`, WT `12.62 ms`.

Capacitance-normalized best-fit rheobase medians are `2.650 pA/pF` in SCA3 and `0.642 pA/pF` in WT.

The working biological interpretation is therefore two-dimensional: SCA3 cells show reduced absolute excitability together with slower suprathreshold spike-train dynamics.

## Robustness to parameter non-identifiability

Parameter non-identifiability does not erase the primary endpoint separation. Endpoint ensemble v1.0.1 contains 64 admissible HR support states from 18 biological cells, with cell-level biological weighting preserved.

The experimental centroid distance in the primary robust-scaled endpoint space `(log10 J_rheo, log10 firing_rate_q75)` is `1.831` pooled-MAD units. Across all enumerated within-cell support-state extremes, the SCA3-WT difference remains positive for normalized rheobase and negative for `q=0.75` firing rate. Hence the primary endpoint direction is robust to the currently identified parameter uncertainty.

Mean ISI at `q=0.75` is retained as a parallel observable because firing rate and ISI are not strictly redundant in finite/transient spike trains.

## Phase-resolved status

Best-fit HR trajectories show phase-localized differences, but full phase profiles are substantially less robust to near-optimal HR alternatives than rheobase, firing rate, and ISI. Absolute `z(phi)` is particularly unstable. Phase-resolved geometry is therefore exploratory/supporting at this stage and must be uncertainty-aware before being used as a primary mechanistic claim.

## Statistical scope

Only two recovered animals are currently available per group, and several WT animal IDs remain unresolved. Formal animal-level genotype p-values are therefore not reported. Cell-level effects describe the recorded-cell ensemble and recovered-animal summaries are descriptive.

## Next frozen boundary

The repository is now transition-ready. The next analysis should build weighted WT→SCA3 trajectory ensembles from biological cell pairs and their admissible within-cell HR support states, while checking staging in both rheobase+firing-rate and rheobase+ISI observable projections.
