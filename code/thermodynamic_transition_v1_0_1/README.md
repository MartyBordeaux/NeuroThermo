# NeuroThermo thermodynamic transition v1.0.1

Corrected stochastic/information-thermodynamic overlay on the frozen uncertainty-aware WT→SCA3 transition geometry. No Hindmarsh–Rose parameters are refit and the geometric stage boundaries are not redefined.

The scientific server profile uses all **264 support-state scenarios** belonging to the **32 core-secure biological WT×SCA3 pairs**. Near-optimal HR solutions are within-pair model uncertainty, not independent biological replicates.

Primary numerical integration uses `dt = 0.05 ms`; a hard completeness gate requires all selected support-state scenarios. Post-processing is streaming and does not materialize millions of trajectory rows.

The stationary NESS layer computes Shannon entropy, endpoint KL divergences, Fisher information, model EPR and density diagnostics. Finite-state Hatano–Sasa exact propagation is the mathematical identity gate. Markov Monte Carlo, continuous trajectories, SMC, path IFT and Crooks histogram symmetry are retained as empirical rare-event diagnostics. Classical Jarzynski/work-Crooks remain disabled because no physical `H(x;lambda)`, `W[path]` or calibrated `beta` is frozen.

Run from this directory with the frozen input tables available under `frozen/`:

```bash
python3 -m thermo_v1_0_1.cli validate --config configs/server_thermodynamic_v1_0_1.yaml
./run_server.sh configs/server_thermodynamic_v1_0_1.yaml
./run_noise_sensitivity.sh
./run_dt_convergence.sh
```

The original pipeline artifact is `neurothermo_thermodynamic_transition_v1.0.1.zip`; SHA256 recorded for this frozen version: `084e0a8345ccf1f7d23b2f3305241cf1584ce054a08375218c08df23e9126080`.
