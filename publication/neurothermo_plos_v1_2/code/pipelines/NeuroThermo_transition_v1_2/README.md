# NeuroThermo transition ensemble v1.2

Two-dimensional WT→SCA3 transition surface over independent intrinsic and drive progress coordinates.

This stage starts from the frozen endpoint ensemble and the corrected v1.1 transition geometry. It does **not** refit any cell.

## Coordinates

`p_intrinsic` controls the intrinsic HR coordinates:

- `b`: linear WT→SCA3 interpolation;
- `s`: linear interpolation;
- `r`: log-linear interpolation.

`p_drive` controls the input-related coordinates:

- `kappa_I`: log-linear interpolation;
- applied q=.75 protocol current `J`: linear interpolation.

The active support window follows `p_intrinsic`, exactly preserving the convention used by transition v1.0 along the coupled path. This is a protocol-analysis choice, not a claim that support duration is an intrinsic HR parameter.

The previous paths are embedded in this plane:

- coupled: `p_drive = p_intrinsic`;
- drive early: `p_drive = 1-(1-p_intrinsic)^2`;
- drive late: `p_drive = p_intrinsic^2`.

## Primary staging geometry

Primary projection: frozen corrected ISI geometry from transition v1.1.

- WT-exit: `A_ISI = 0.135829`;
- balance: `A_ISI = 0.5`;
- SCA3-entry: `A_ISI = 0.797856`.

The corrected experimental active-rate projection is retained as an independent secondary map.

## Biological/model weighting

The full server profile contains 72 WT×SCA3 biological cell pairs and 988 admissible parameter-support scenarios. Biological pairs are equally weighted; near-optimal parameter solutions divide the weight within each pair and are not treated as additional biological replicates.

## Server run

```bash
unzip neurothermo_transition_ensemble_v1.2.0.zip
cd neurothermo_transition_ensemble_v1_2

python3 -m transition_v1_2.cli validate \
  --config configs/server_transition_v1_2.yaml

./run_server.sh configs/server_transition_v1_2.yaml
```

No editable install is required when run from the unpacked project directory.

The full profile is a 31×31 surface for all 988 support scenarios: 949,468 simulated states. Checkpoints are written per support scenario and `run_server.sh` resumes them automatically.

After the run, package the complete result directory without losing its top-level folder:

```bash
./pack_results.sh
```

This creates `results_transition_ensemble_v1_2.zip` containing the single root folder `results_transition_ensemble_v1_2/`.
