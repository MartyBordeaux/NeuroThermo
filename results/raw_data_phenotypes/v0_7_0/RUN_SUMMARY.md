# NeuroThermo v0.7.0 — predictive-dynamics mechanistic validation

## Analysis unit and cohort

Each cell is an independent analysis unit. Animal identifiers are removed before analysis. The frozen cohort contains 13 WT and 7 SCA3 cells on the 100–600 pA current grid.

## Primary adjustment

Log predictive information is predicted separately at each current by label-blind leave-one-cell-out ridge regression. The primary covariates are log1p firing rate, log mean ISI with an explicit missingness indicator, log baseline noise, and stationary-sample count. WT/SCA3 labels are never used by the residualizer. The primary ridge penalty is lambda=1.0 and was frozen before exact group tests.

## Primary results

- Adjusted two-domain burden curve AUC difference: 1.28058; exact p=0.000386997.
- Adjusted current-wise burden maxT exact p: 0.029773.
- Residual predictive-dynamics curve AUC difference: 1.28658; exact p=0.00216718.
- Adjusted exits: WT 0/13; SCA3 2/7.
- Exiting cells: SCA3_02=500 pA, SCA3_09=500 pA.
- Adjusted I_exit restricted-mean difference WT minus SCA3: 42.8571 pA; exact p=0.00122549.

## Ridge sensitivity

lambda 0.1: WT 0, SCA3 0 (none); lambda 1.0: WT 0, SCA3 2 (SCA3_02;SCA3_09); lambda 10.0: WT 0, SCA3 1 (SCA3_02).

## Interpretation boundary

Persistence after adjustment supports a predictive-dynamics contribution not explained by the included activity and technical covariates. It does not establish causal mechanism, statistical independence of domains, disease time, irreversibility, or a thermodynamic phase transition. I_exit remains a threshold under imposed current stress.
