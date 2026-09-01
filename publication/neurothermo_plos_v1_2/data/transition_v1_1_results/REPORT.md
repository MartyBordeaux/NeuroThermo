# Transition ensemble v1.1 — correction/reprojection report

## Scope

v1.1 reuses all 121,524 HR transition states from v1.0 and performs **zero new HR simulations**. The only correction is semantic/geometric: v1.0 `support_rate_hz` is an active-window rate measured in a fixed window starting at the first model spike, so v1.1 names it `active_support_rate_hz` and no longer projects it against onset-normalized firing rate.

The primary staging projection remains the latency-invariant ISI projection from v1.0 and is reproduced numerically to machine precision (maximum absolute difference in A: 8.038e-14). The secondary projection uses experimental `active_rate_hz` at q=0.75 from dynamic v2.1.

## Corrected endpoint boundaries

Primary ISI projection:

- WT-exit: A = 0.135829;
- balance: A = 0.5;
- SCA3-entry: A = 0.797856.

Corrected active-rate projection:

- WT-exit: A = 0.148329;
- balance: A = 0.5;
- SCA3-entry: A = 0.683282.

The active-rate reference clouds do not overlap.

## Primary ISI staging

Values are medians [Q25-Q75] across the 32 core-secure biological WT×SCA3 pairs after collapsing parameter-support uncertainty within each pair.

| path family | WT-exit | balance | SCA3-entry |
|---|---:|---:|---:|
| drive early | 0.234 [0.148-0.356] | 0.523 [0.403-0.599] | 0.717 [0.617-0.830] |
| coupled | 0.398 [0.105-0.490] | 0.676 [0.590-0.707] | 0.837 [0.773-0.869] |
| drive late | 0.401 [0.081-0.656] | 0.793 [0.710-0.830] | 0.878 [0.743-0.930] |

Thus the principal result from v1.0 survives v1.1 unchanged: advancing the input-related drive moves the macroscopic transition earlier, whereas delaying it moves balance and SCA3-entry later.

## Independent active-rate check

Corrected active-rate stage medians in the same core-secure subset are:

| path family | WT-exit | balance | SCA3-entry |
|---|---:|---:|---:|
| drive early | 0.276 [0.159-0.365] | 0.476 [0.352-0.585] | 0.606 [0.427-0.655] |
| coupled | 0.409 [0.167-0.521] | 0.683 [0.458-0.753] | 0.774 [0.582-0.836] |
| drive late | 0.426 [0.087-0.690] | 0.775 [0.395-0.850] | 0.820 [0.323-0.908] |

The strongest agreement is at **balance**: the median active-vs-ISI stage-position gap is about 0.02 in every path family, and the within-pair support-weight agreement at the adopted 0.10 tolerance has median 1.0 for all three families. WT-exit is also close. SCA3-entry is less interchangeable, particularly for drive-early, so v1.1 does not use an active/ISI average as the primary stage marker.

## Occupancy

At the ISI-defined stages, median model occupancy in core-secure pairs remains approximately 0.98-0.99 across path families. Thus the primary staging result is not generated mainly by isolated one- or two-spike events within the active support window.

## Audit

All correction checks pass. The best endpoint active-window rate agrees with independently frozen dynamic-v2.1 model `active_rate_hz` with median symmetric relative difference 1.13%.

## Interpretation

v1.1 validates the ISI-based staging used in v1.0 and repairs the secondary rate geometry. The next analysis can therefore treat ISI staging as frozen geometry and use corrected active rate as an independent check before constructing the two-dimensional intrinsic-progress × drive-progress sensitivity map.
