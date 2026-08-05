# Server run protocol

Raw ABF directories remain unchanged:

```text
~/neurothermo/WT
~/neurothermo/SCA3
```

This pipeline consumes frozen derived endpoints, not ABF files. The frozen
table is excluded from Git and remains in protected server storage. Before
installing the package, stage it locally:

```bash
mkdir -p src/neurothermo_per_cell/data
install -m 600 "$HOME/neurothermo_data/frozen_v2_w20_observations.csv" \
  src/neurothermo_per_cell/data/frozen_v2_w20_observations.csv
```

The GitHub Actions workflow uses this same default source path; set the
repository variable `NEUROTHERMO_BENCHMARK` only when the protected file is
stored elsewhere on the self-hosted runner.

## Smoke test

```bash
cd "$HOME/neurothermo/per_cell_pipeline"
source .venv/bin/activate
PYTHONPATH=src python -m neurothermo_per_cell.cli run --config configs/smoke.yaml
cat runs/smoke/run_metadata.json
cat runs/smoke/failures.csv
```

## Preliminary run

```bash
nohup bash -lc '
cd "$HOME/neurothermo/per_cell_pipeline"
source .venv/bin/activate
PYTHONPATH=src python -m neurothermo_per_cell.cli run --config configs/preliminary.yaml
' > "$HOME/neurothermo/per_cell_preliminary.log" 2>&1 &
echo $! > "$HOME/neurothermo/per_cell_preliminary.pid"
```

Monitor:

```bash
tail -n 30 "$HOME/neurothermo/per_cell_preliminary.log"
find runs/preliminary/fits -name '*.json' | wc -l
```

Expected: 192 fit JSON files and 385 lines in `result_rows.csv`, including the
header. Inspect:

```text
runs/preliminary/run_metadata.json
runs/preliminary/failures.csv
runs/preliminary/summary/full_cell_fits.csv
runs/preliminary/summary/within_cell_cv_metrics.csv
runs/preliminary/summary/cell_parameters_long.csv
runs/preliminary/summary/convergence_summary.csv
```

The configuration uses `workers: 4`. Set it to the number of CPU cores allocated
by the scheduler, not the total number physically present on a shared server.

## Full run

The supplied full profile has 960 tasks, five seeds, `dt=0.1 ms`, and a larger
optimizer budget. Start it only after a stop/go review of preliminary results.

```bash
nohup bash -lc '
cd "$HOME/neurothermo/per_cell_pipeline"
source .venv/bin/activate
PYTHONPATH=src python -m neurothermo_per_cell.cli run --config configs/full.yaml
' > "$HOME/neurothermo/per_cell_full.log" 2>&1 &
echo $! > "$HOME/neurothermo/per_cell_full.pid"
```
