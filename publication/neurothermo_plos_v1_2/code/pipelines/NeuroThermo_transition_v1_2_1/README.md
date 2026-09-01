# NeuroThermo transition ensemble v1.2.1

Scenario-first uncertainty correction for the frozen v1.2 intrinsic × drive surfaces.

This stage performs **zero new HR simulations**. It reads the 988 existing `scenario_XXXX.csv.gz` checkpoint surfaces from `results_transition_ensemble_v1_2` and corrects the order of uncertainty aggregation:

1. detect each stage crossing separately for every admissible WT→SCA3 support-state scenario;
2. compute v1.1-compatible weighted crossing quantiles over the support states in which a marker exists, while reporting the total crossing-support weight separately;
3. aggregate the marker distribution within each biological pair;
4. aggregate biological pairs with equal pair weight;
5. compute drive sensitivity per scenario before any within-pair or ensemble aggregation.

The primary projection remains the frozen v1.1 ISI geometry. Active-rate remains secondary.

## Server run

Keep either `results_transition_ensemble_v1_2/` or `results_transition_ensemble_v1_2.zip` next to this pipeline folder.

```bash
python3 -m transition_v1_2_1.cli validate --config configs/server_transition_v1_2_1.yaml
./run_server.sh configs/server_transition_v1_2_1.yaml
```

No `pip install -e .` is required when commands are run from this folder.

After completion:

```bash
./pack_results.sh
```
