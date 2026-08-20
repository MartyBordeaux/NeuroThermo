# Validation report: v3.9.0

Before scientific use, v3.9 must satisfy:

1. Same frozen cohort and spike events as v3.8.
2. Exactly one scientific bound change: `s_min = 0.05`; `s_max = 15`, and all `b`, `r`, `kappa_I` bounds unchanged.
3. Same bounds for WT and SCA3.
4. Exact first-spike additive alignment; no time rescaling or last-spike anchor.
5. Alignment never changes model spike count.
6. Same binary rheobase constraint and no plateau/ISI/adaptation loss.
7. Full global search for every cell, with v3.8 solutions only as seeds.
8. Identifiability separated alternatives use the original v3.6 reference separation, not a fraction of the enlarged v3.9 domain.
9. `s_boundary_stress_summary.csv` must report whether each optimum falls below the old 0.25 bound and proximity to the new 0.05 bound.
10. SCA3_05 receives no special rescue, exclusion, or sensitivity branch.
