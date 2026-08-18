# Transition analysis freeze: v1.1 → v1.3.1

## v1.1 — corrected staging geometry

v1.1 reprojects the v1.0 transition states without new HR simulations. The primary observable is latency-invariant mean ISI combined with model rheobase. The three frozen ISI thresholds are:

- WT-exit: `A=0.1358293233470019` (maximum A among core-secure WT endpoints);
- balance: `A=0.5` (midpoint of the WT→SCA3 centroid axis);
- SCA3-entry: `A=0.7978563373093712` (minimum A among core-secure SCA3 endpoints).

For 32 core-secure biological WT×SCA3 pairs, the coupled-path medians are WT-exit `0.3983`, balance `0.6764`, SCA3-entry `0.8370`. Early drive shifts the stages earlier; late drive delays balance and SCA3-entry.

The secondary active-rate projection independently supports the transition ordering but is not averaged with ISI into the primary marker.

## v1.2.1 — scenario-first intrinsic × drive surface

The 2D map uses independent coordinates for intrinsic progress `(b,r,s)` and combined drive progress `((kappa_I), J)`. v1.2.1 fixes the uncertainty-ordering issue in v1.2 by finding each support-scenario crossing first, then aggregating within a biological pair, then across biological pairs.

In the core-secure subset, median drive-dominance fraction rises from `0.339` at WT-exit to `0.403` at balance and `0.518` at SCA3-entry. Thus drive is not uniformly dominant; its relative contribution increases toward the SCA3-like region.

The old one-dimensional paths are recovered closely from the 2D surface except for `drive_late WT-exit` (`|Δp|≈0.112`), which is retained as a 31×31 grid-resolution diagnostic. The frozen v1.1 one-dimensional estimate remains authoritative for that particular crossing.

## v1.3 — factorized drive decomposition

The combined drive is decomposed into:

- `kappa_I`: the fitted HR input-scaling coordinate;
- `J`: the experimenter-imposed current protocol;
- a non-additive interaction term.

At core-secure ISI stage bands, the median absolute `kappa_I` effect exceeds the `J` effect at WT-exit, balance, and SCA3-entry. Median signed interactions are modest and negative (`-0.040`, `-0.025`, `-0.022`, respectively), supporting weak sub-additivity rather than strong positive synergy.

The main mechanistic interpretation is therefore not that external current is a disease driver. Instead, the fitted input-scaling coordinate contributes strongly to SCA3-directed movement while the experimental `J` protocol partially offsets that movement. `J` remains protocol sensitivity only. Because raw `kappa_I` covaries strongly with capacitance, it is not interpreted as a direct membrane conductance or channel-level mechanism.

## v1.3.1 figure correction

Numerical v1.3 results are unchanged. The original three boundary panels were blank because the plotting function used `b.mode`, which pandas resolves to `DataFrame.mode()` instead of the `mode` column. v1.3.1 changes the filter to explicit column access (`b['mode']`, `b['stage']`) and regenerates the plots.

## Frozen interpretation

The WT→SCA3 transition is staged and uncertainty-aware. Intrinsic and input-related coordinates jointly determine the transition. Relative input sensitivity increases toward late SCA3-like progression, but the factorized analysis attributes most of this late drive effect to the model input-scaling coordinate rather than the imposed current protocol.

The transition analysis is now ready for a thermodynamic overlay rather than another geometric redefinition.
