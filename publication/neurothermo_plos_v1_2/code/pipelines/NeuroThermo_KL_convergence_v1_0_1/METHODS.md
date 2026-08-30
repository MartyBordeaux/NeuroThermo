# Prespecified KL-convergence analysis

The full frozen core-secure ensemble is simulated at integration steps 0.05,
0.025, and 0.0125 ms with the five previously frozen seeds. Burn-in, sampling
duration, sampling stride, noise amplitudes, path positions, density estimator,
and endpoint-relative KL definition are otherwise unchanged.

For a scenario and path position, all three integration steps consume the same
finest-grid Gaussian stream. Coarser Brownian increments are sums of adjacent
finest-grid increments. This common-random-number construction separates
integration-step effects from avoidable Monte Carlo differences.

The calculation uses two deterministic passes. In the pilot pass, every one of
the 15 step--seed paths is simulated and only the coordinate-wise sample
extrema are retained. For each scenario, a single reference histogram grid is
then constructed from the combined extrema with a 2% margin. In the analysis
pass, the same paths are regenerated from the same seeds and evaluated on this
grid. Any task retaining less than 100% of its samples terminates the run. Raw
trajectories are not stored.

For each view, the endpoint-affinity curve is

\[
\Delta D_{KL}(p)=D_{KL}[\rho_p\Vert\rho_{WT}]
                 -D_{KL}[\rho_p\Vert\rho_{SCA3}].
\]

Zero therefore denotes equal endpoint-relative divergence. Six frozen marker
variants are evaluated: the isotonic, first, and persistent crossings of the
seed-median curve and the 25th, 50th, and 75th percentiles of seed-specific
isotonic crossings.

Two noncommuting analysis orders are compared:

1. marker-first: extract each support-scenario marker and then take the
   within-pair weighted median;
2. curve-first: aggregate support-scenario KL curves within each pair and then
   extract the marker.

Pair results are subsequently summarized by endpoint cell and by crossed WT
animal--SCA3 animal-day stratum. Leave-one-animal-out analyses omit each of the
four WT animals and two SCA3 animal-days in turn.

The formal verdict is hierarchical. Fatal gates cover task completeness, full
grid coverage, `xyz` scenario convergence, the fine-step pair- and
cell-balanced effect for the seed-median isotonic marker, step-size stability,
and agreement between aggregation orders. Failure of any fatal gate yields
`REMOVE_KL_RESULT`. The `xy` view, alternative marker variants, animal-pair
summaries, and leave-one-animal-out checks are supporting because they test
representation and dependence sensitivity rather than the primary full-state
claim. Failure confined to these gates yields
`KEEP_AS_ENSEMBLE_RESULT_WITH_LIMITATIONS`; passage of every gate yields
`KEEP_AS_MAIN_RESULT`. Thresholds are frozen in
`configs/server_kl_convergence_v1_0_1.yaml`.
