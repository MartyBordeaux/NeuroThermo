# Validation and decision gates

The release smoke-test executes one scenario, five path positions, one seed,
fixed-point analysis, current continuation, density/current estimation, Fisher
geometry, Markov analysis, both fluctuation diagnostics, output assembly, and
figure creation.

The full configuration additionally requires exactly 264 support scenarios,
32 dependent cell-pair combinations, four WT animals, two SCA3 animals, and
eight animal-pair combinations. Input discovery is forbidden: the operator must
pass the exact frozen directory.

The formalism gate requires at least 95% of scenario-seed-path Markov estimates
to have detailed-balance violation <= 0.05. Failure produces a `NESS` verdict.

The continuous-current branch is considered numerically adequate only when its
relative stationarity-divergence residual is <= 0.25. An inadequate branch is
labelled `DIAGNOSTIC_INVALID` and cannot support either equilibrium or NESS. It
does not override a decisive independent Markov time-reversal failure.

These thresholds are intentionally strict screening thresholds, not p-values.
The sensitivity table must also be inspected before interpreting a borderline
result. Markov lags of 2.5, 10, and 20 ms remain required sensitivity outputs.

Protocol validation requires every schedule to have exactly
`protocol.n_points` unique positions, including both endpoints. Linear and
adaptive schedules must therefore have identical state counts.

Cycle-affinity rows must include `scenario_id`, `biological_pair_key`, `seed`,
and `p`. Missing provenance is a packaging failure.

The classical Jarzynski–Crooks gate is conjunctive. Passing detailed balance is
insufficient; a physical beta, energy, and conjugate work mapping for all five
morphed coordinates are also required. The shipped configuration keeps this
mapping disabled.
