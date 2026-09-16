# NeuroThermo transition ensemble v1.3

Factorized decomposition of the frozen v1.2 combined drive into:

- fitted HR input scaling `kappa_I`;
- experimental current protocol `J`.

The analysis does **not** refit any cell. It reuses the complete combined `(intrinsic, drive)` surface from transition v1.2 and computes two additional 2D surfaces for every selected WT→SCA3 support-state scenario.

## Factorial design

At fixed intrinsic progress `p_i` and component progress `q`:

- **combined**: `kappa_I=q`, `J=q` (reused from v1.2);
- **kappa_only**: `kappa_I=q`, `J=p_i`;
- **J_only**: `kappa_I=p_i`, `J=q`;
- **coupled baseline**: `kappa_I=p_i`, `J=p_i`.

The intrinsic coordinates `(b,r,s)` always follow `p_i`.

For either projected observable `A`, the decomposition is

`K = A(kappa=q,J=p_i) - A_base`

`J_effect = A(kappa=p_i,J=q) - A_base`

`Combined = A(kappa=q,J=q) - A_base`

`Interaction = Combined - K - J_effect`.

The interaction therefore measures non-additivity of the two drive components at matched progress; it is not a Bayesian interaction coefficient.

## Important interpretation guardrails

`J` is experimentally imposed current density. The J-only surface is therefore **protocol sensitivity**, not a disease trajectory.

`kappa_I` is a fitted HR input-scaling coordinate. Its raw value strongly covaries with capacitance in the frozen cohort, so the kappa-only surface is a model-coordinate sensitivity analysis and must not be presented as an independently established biological SCA3 mechanism.

The primary geometry remains the frozen v1.1 ISI projection; corrected active-rate is secondary validation.

## Server run

Place either `results_transition_ensemble_v1_2/` or `results_transition_ensemble_v1_2.zip` next to this project directory.

```bash
cd ~/neurothermo/v4
unzip neurothermo_transition_ensemble_v1.3.0.zip
cd neurothermo_transition_ensemble_v1_3

python3 -m transition_v1_3.cli validate \
  --config configs/server_transition_v1_3.yaml

./run_server.sh configs/server_transition_v1_3.yaml
```

The server profile uses all 988 support-state scenarios and a `31 x 31` grid. It adds two new surfaces, so it evaluates `1,898,936` new HR states. The original combined surface is not recomputed.

The run is resumable at the scenario/mode checkpoint level.

After completion:

```bash
./pack_results.sh
```

This creates `results_transition_ensemble_v1_3.zip` containing one top-level results directory.

## Preliminary profile

```bash
./run_server.sh configs/preliminary_transition_v1_3.yaml
```

This uses only the 72 best×best biological scenarios on a `21 x 21` grid and is intended only for workflow/performance validation.
