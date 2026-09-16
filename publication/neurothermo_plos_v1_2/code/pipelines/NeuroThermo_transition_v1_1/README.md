# NeuroThermo WT→SCA3 transition ensemble v1.1

This is a **reprojection/correction stage** for the already completed transition ensemble v1.0. It performs **zero new Hindmarsh–Rose simulations**.

## Why v1.1 exists

In v1.0 the synthetic transition quantity `support_rate_hz` was measured in a fixed window starting at the **first model spike**, so it is a latency-invariant active-window rate. However, the v1.0 rate projection was referenced to the older onset-normalized `firing_rate_hz`. These are not the same observable.

v1.1 therefore:

- renames the semantic quantity to `active_support_rate_hz` without changing a single simulated value;
- retains the v1.0 **ISI projection unchanged as the primary staging projection**, because mean ISI was already latency-invariant and internally consistent;
- builds a new independent active-rate reference from the experiment-supported `q=0.75` `active_rate_hz` values frozen from dynamic characterization v2.1;
- recomputes active-rate projection coordinates and stage markers;
- reports active-vs-ISI agreement instead of treating their average as the primary result;
- summarizes occupancy fraction at the ISI-defined stages as an additional latency-invariant descriptor.

## Stage definitions

The primary ISI staging thresholds are frozen from v1.0:

- WT-exit: maximum `A_ISI` among core-secure WT endpoints;
- balance: `A_ISI=0.5`, the midpoint between WT and SCA3 centroids;
- SCA3-entry: minimum `A_ISI` among core-secure SCA3 endpoints.

The active-rate projection uses the same definitions but a newly fitted reference cloud from experimental `active_rate_hz` at `q=0.75`.

No cell is excluded from the full uncertainty ensemble. `core_secure_pairs` remain a stringent reference subset for stage summaries.

## Run

Place this folder next to either the unzipped `results_transition_ensemble_v1_0/` directory or `results_transition_ensemble_v1_0.zip`.

```bash
python3 -m transition_v1_1.cli validate --config configs/server_transition_v1_1.yaml
./run_server.sh configs/server_transition_v1_1.yaml
```

No editable install is required.
