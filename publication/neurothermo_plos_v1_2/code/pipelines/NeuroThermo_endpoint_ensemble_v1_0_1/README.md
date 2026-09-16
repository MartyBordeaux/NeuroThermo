# NeuroThermo endpoint ensemble v1.0

This stage freezes uncertainty-aware WT and SCA3 endpoint ensembles for the next WT→SCA3 transition analysis. It does **not** refit Hindmarsh–Rose parameters and does not perform a transition calculation.

## Scientific definition

All 18 accepted multi-sweep cells are retained (12 WT, 6 SCA3). Every cell contributes one best v3.9 HR parameter vector and, where identifiability analysis found them, the near-optimal alternative vectors. Alternatives are treated as a **discrete model-uncertainty support**, not as posterior samples and not as additional biological replicates.

To prevent non-identifiable cells from being over-weighted, each biological cell has equal weight within its group; that weight is divided equally among the available HR support members for that cell. These are numerical support weights, not probabilities.

The complete observable record is

`(J_rheo, firing_rate(q=.75), mean_ISI(q=.75), firing_rate(q=.50), mean_ISI(q=.50))`.

No missing q=.50 value is imputed. q=.75 is experimentally supported for all 18 cells; q=.50 is supported for 17/18 cells.

For transition geometry the primary common anchor is intentionally reduced to two non-redundant coordinates:

`log10(J_rheo)` and `log10(firing_rate(q=.75))`.

Mean ISI is retained as an independent validation/readout but is not counted as a second rate-like geometric axis, because rate and ISI are strongly reciprocal. The two core coordinates are robustly standardized by the pooled best-cell median and MAD. Actual HR parameters remain attached to every endpoint support member so that later transition paths can be generated in model space while remaining anchored to an observable phenotype space.

## Outputs

- `endpoint_cells_full_observable.csv`: full experimental/model endpoint record and dynamic-security flags for every primary cell.
- `endpoint_cells_transition_core.csv`: core transition coordinates for the 18 best-cell endpoints.
- `endpoint_solution_support.csv`: best + near-optimal HR support members with observable predictions.
- `transition_ready_endpoint_support.csv`: compact table intended as the input to the transition pipeline.
- `transition_core_transform.csv`: frozen robust coordinate transformation.
- `group_endpoint_summary.csv`: descriptive group medians/IQR and Cliff effects; no animal-level p-values.
- `uncertainty_decomposition.csv`: between-cell biological spread versus within-cell parameter-solution spread.
- `endpoint_geometry.csv`: WT and SCA3 core centroids and their descriptive distance.

## Reproducibility

```bash
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -e .
python3 -m endpoint_v1.cli validate --config configs/server_endpoint_v1_0.yaml
./run_server.sh configs/server_endpoint_v1_0.yaml
```

The animal-ID layer is preserved as metadata. Because only two recovered animals are available per genotype and several WT cells remain unresolved, this stage performs no genotype-level animal p-values.
