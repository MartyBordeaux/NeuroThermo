# Methods and definitions

## Frozen stochastic model

For state `X=(x,y,z)` and diagonal diffusion `D`, the pipeline integrates

```text
dx = [y - a x^3 + b x^2 - z + kappa_I J] dt + sqrt(2 D_x) dW_x
dy = [c - d x^2 - y] dt                         + sqrt(2 D_y) dW_y
dz = [r(s(x-x_R)-z)] dt                         + sqrt(2 D_z) dW_z
```

The constructed coordinate `p` interpolates `b,s,J` linearly and `r,kappa_I`
log-linearly between each frozen WT/SCA3 endpoint. It is model-space position,
not disease time.

## Potential, current, and detailed balance

The empirical nonequilibrium potential is

`phi_p(X) = -log rho_ss(X|p)`.

For additive diffusion the stationary probability current is

`j_p = F_p rho_p - D grad(rho_p)`,

or `j_p/rho_p = F_p - D grad(log rho_p)`. The pipeline reports its squared
`D^{-1}` norm (an EPR-like proxy), its fraction relative to gradient score
energy, and the relative divergence residual. Histogram estimates are repeated
across two bin/smoothing variants.

Independently, trajectories are mapped to a fixed quantile partition of state
space. At each `p`, a lagged transition matrix is estimated and its invariant
distribution is solved numerically. Detailed-balance violation, entropy per lag,
and supported triangular cycle affinities are reported at the primary lag and
three alternative lags.

These are numerical diagnostics. They do not attach physical energy units or a
temperature to the phenomenological HR state variables.

The continuous-current estimator is declared interpretable only when the
relative divergence residual is at most 0.25. Otherwise it remains a QC output
and is excluded from the equilibrium/NESS decision, which is based on the
independent coarse Markov time-reversal test.

## Fisher metrics

The path Fisher metric is calculated as the centered score variance

`g(p) = Var_rho[partial_p log rho_p(X)]`.

It is calculated for the full `xyz` density, the fast `xy` marginal, and the `z`
marginal. `g_xyz-g_xy` is reported as the conditional slow-variable contribution
from the Fisher chain rule. The local check compares adjacent-state KL divergence
with `0.5 g dp^2`.

The state-space score energy is a distinct quantity:

`I_state,D = integral rho_p [grad log rho_p]^T D [grad log rho_p] dX`.

It diagnoses spatial sharpness and must not be substituted for `g(p)`.

## Thermodynamic length and friction

Path length is `L = integral sqrt(g(p)) dp`. The generalized Hatano–Sasa force is
`X_p(t)=partial_p phi_p(X_t)`. Its variance equals the path FI in the continuum
limit. Version 1.0.1 estimates the normalized positive-lobe correlation time
from the sampled force series and multiplies it by the centered density-based
`g(p)`. This makes covariance amplitude and path FI internally consistent. The
raw sample-to-grid variance ratio remains available as QC.

All compared schedules contain 15 strictly increasing positions drawn from the
same 31-position path. A dynamic-programming selector minimizes mismatch to
equal thermodynamic-length targets while forcing both endpoints and preventing
repeated nearest-neighbour states.

## Fluctuation relations

For a discrete protocol the excess functional is

`Y = sum_k [phi_{p_{k+1}}(X_k)-phi_{p_k}(X_k)]`.

The exact finite-state tilted propagation checks `<exp(-Y)>=1`. Monte Carlo
estimates, standard errors, and exponential-weight ESS are reported separately.
A linear schedule and schedules discretized from the `xyz`-FI, `xy`-FI, and
correlation-aware friction lengths are compared directly.
A second diagnostic uses the log ratio of explicitly defined forward and reverse
coarse-grained path probabilities. Neither quantity is called physical work.

## Biological aggregation

Support scenarios are first averaged with frozen within-cell-pair weights.
Cell-pair summaries are then averaged within each recovered animal pair, and
animal-pair distributions are summarized without treating cells from the same
animal as independent. In the core-secure cohort this gives four WT animals,
two SCA3 animal-days, and eight crossed animal-pair combinations. These are
sensitivity summaries; two SCA3 animals remain insufficient for population-
level genotype inference.

Classical Jarzynski/Crooks is blocked unless all of the following are supplied
and pass: detailed balance, physical inverse temperature, a physical energy
definition, and conjugate work for **every** morphed control (`b,r,s,kappa_I,J`).
Counting only `integral x dJ` while morphing intrinsic parameters is explicitly
rejected by configuration validation.
