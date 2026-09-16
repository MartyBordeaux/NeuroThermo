# NeuroThermo WT→SCA3 transition ensemble v1.0

This pipeline is the first transition/staging calculation built on the frozen NeuroThermo endpoint ensemble v1.0.1. It does **not** refit any cell.

## Scientific unit of analysis

The primary biological layer contains 18 multi-sweep cells (12 WT, 6 SCA3), yielding 72 WT×SCA3 biological cell pairs. The endpoint uncertainty layer contains 64 admissible HR support states (38 WT, 26 SCA3), yielding 988 support-state pair scenarios.

Near-optimal HR solutions are model uncertainty, not additional biological replicates. Each biological pair has equal weight; its weight is divided across the Cartesian product of the two cells' admissible support states.

## Transition path

For a WT support state and an SCA3 support state, path progress is `p∈[0,1]`.

- `b` and `s` are interpolated linearly.
- `r` is interpolated linearly in `log(r)`.
- `kappa_I` is interpolated linearly in `log(kappa_I)`.
- the applied current is anchored to the experimentally supported q=0.75 current of each endpoint cell.
- the active observation window is interpolated between endpoint q=0.75 experimental-support windows.

Three path families are implemented:

1. `coupled`: intrinsic parameters and drive progress together;
2. `drive_early`: `kappa_I` and applied current progress earlier than `(b,r,s)`;
3. `drive_late`: `kappa_I` and applied current progress later than `(b,r,s)`.

The last two are perturbation/protocol-sensitivity paths, not claims about the biological order of degeneration.

## Why latency is not reintroduced

At each synthetic intermediate state, absolute excitability is characterized separately by the model rheobase. Suprathreshold dynamics are measured in a window starting at the **first model spike**, so onset latency remains a nuisance coordinate, consistent with the v3.6-v3.9 fit logic. The active window duration comes from the experiment-supported q=0.75 region of the two endpoint cells.

The transition rate metric is therefore named `support_rate_hz`; it is not silently equated with the original raw firing-rate definition. Endpoint protocol validation is written to `endpoint_protocol_validation.csv`.

## Staging geometry

Two projections are retained in parallel:

- rheobase + supported rate;
- rheobase + mean ISI.

Reference centroids and robust scaling are estimated from the best-fit, `core_q75_secure` endpoint cells. For each path point the pipeline calculates the directed WT→SCA3 coordinate `A` and orthogonal deviation from the endpoint axis.

Stage markers are:

- `WT-exit`: first persistent crossing of the upper boundary of the stringent WT reference envelope (`max A` by default);
- `balance`: first persistent crossing of `A=0.5`;
- `SCA3-entry`: first persistent crossing of the lower boundary of the stringent SCA3 reference envelope (`min A` by default).

Markers are reported independently for the rate and ISI projections. A consensus marker is created only when the two projections agree within the configured tolerance.

`linear`, `smoothstep`, `quadratic`, and `sqrt` schedule morphs are reported by mapping the geometric path-progress marker `p` back to a schedule coordinate `u`. These morphs are reparameterizations of a given path, not different geometric paths.

## Server run

From the extracted project directory:

```bash
python3 -m transition_v1.cli validate --config configs/server_transition_v1_0.yaml
./run_server.sh configs/server_transition_v1_0.yaml
```

The full server configuration evaluates all 988 support scenarios, three path families, and 41 points per path. Checkpoints are written per scenario×family, so interrupted runs can be resumed by rerunning the same command.

A quicker best-fit-only run is available:

```bash
./run_server.sh configs/preliminary_transition_v1_0.yaml
```

## Interpretation boundary

This v1.0 stage tests whether uncertainty-aware HR paths show reproducible geometric staging between the frozen WT and SCA3 endpoint clouds. Fisher information, KL balance, entropy, EPR, Hatano-Sasa, or other thermodynamic quantities are deliberately **not** included yet. They should only be added after geometric staging survives endpoint uncertainty and path-family sensitivity.
