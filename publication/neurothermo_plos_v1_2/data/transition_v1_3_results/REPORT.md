# Transition ensemble v1.3 — factorized drive decomposition

Combined drive is decomposed into fitted HR input scaling kappa_I and experimental current protocol J. The frozen v1.2 combined surface is reused; only kappa-only and J-only maps are newly simulated.

Factorial contrasts at fixed intrinsic progress p_i and matched component progress q:

- K = A(kappa=q, J=p_i) - A(kappa=p_i, J=p_i)
- J = A(kappa=p_i, J=q) - A(kappa=p_i, J=p_i)
- Combined = A(kappa=q, J=q) - A(kappa=p_i, J=p_i)
- Interaction = Combined - K - J

The J-only surface is protocol sensitivity, not a disease parameter trajectory. Raw kappa_I is a fitted model coordinate and remains confounded with capacitance; it is not interpreted as an independent biological phenotype.

## Core-secure ISI interaction near frozen stage boundaries

- SCA3_entry: |K|=0.2781, |J|=0.1513, |interaction|=0.0328, signed interaction=-0.0222.
- WT_exit: |K|=0.3087, |J|=0.2149, |interaction|=0.0322, signed interaction=-0.0396.
- balance: |K|=0.2952, |J|=0.1524, |interaction|=0.0314, signed interaction=-0.0251.
