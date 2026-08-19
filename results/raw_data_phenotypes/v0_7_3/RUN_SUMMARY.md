# NeuroThermo v0.7.3 — rank-Gaussian exact-spectrum Fourier validation

## Frozen design

All 13 WT and 7 SCA3 cells are retained. The cell is the independent unit; animal-level inference is absent. The primary endpoint remains the imported v0.7.2 shuffle-centered non-overlapping PI at 4 ms. Failed v0.7.2 IAAFT results are excluded from inference.

## Primary endpoint

WT-minus-SCA3 cell-AUC difference: 0.14784769 nats; exact expected-direction p=0.0091847265; maxT p across lags=0.0093653251.

## Exact-spectrum sensitivity

The voltage trace is transformed monotonically to empirical normal scores. This leaves every ordinal code unchanged. Fourier phases are randomized while the complete rFFT magnitude spectrum of the rank-Gaussian trace is preserved.

Fourier-centered 4-ms WT-minus-SCA3 cell-AUC difference: 0.1365028 nats; exact expected-direction p=0.0095975232; maxT p across lags=0.0095975232.

All exact-spectrum fidelity gates passed: True. Maximum ordinal-PI difference introduced by rank Gaussianization: 2.22e-16 nats.

## Interpretation boundary

Survival against this null identifies ordinal temporal structure not reproduced by a linear Gaussian-copula process with the same rank-Gaussianized spectrum. It does not prove entropy production, thermodynamic irreversibility, causal mechanism, disease time or a phase transition.
