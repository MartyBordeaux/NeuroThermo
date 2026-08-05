# NeuroThermo per-cell pipeline v0.2.0

This branch replaces the phenotype-level compromise fit with joint fitting of all
current steps from one cell. Parameters are constant across the cell's sweeps.
It does **not** fit individual sweeps separately.

## Scientific outputs

1. `full` fits estimate one parameter vector per cell from all 13 current levels.
2. `cv0`–`cv2` fits hold out interleaved positive-current levels and quantify
   within-cell current-response prediction.
3. Differential-evolution history and the final near-optimal population diagnose
   convergence and local practical non-identifiability.
4. Profile likelihood is a separate post-fit calculation and accepts only
   `CONVERGED` full-cell fits.

The 20-ms observation table is reconstructed exactly from frozen benchmark
v2.0.0 outputs: 312 rows, 15 WT cells, 9 SCA3 cells, 13 current levels per cell.
It contains derived endpoints, not ABF recordings. It is deliberately excluded
from Git. Before local execution, place it at
`src/neurothermo_per_cell/data/frozen_v2_w20_observations.csv`; for a
self-hosted runner, retain it in protected server storage and stage it at that
path before running the workflow. A benchmark CSV can instead be supplied
through `data_path`.

## Installation on Python 3.9.25

```bash
cd "$HOME/neurothermo/per_cell_pipeline"
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade "pip<26" "setuptools<76" wheel
python -m pip install -e .
python -m pip check
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Commands

```bash
PYTHONPATH=src python -m neurothermo_per_cell.cli validate --config configs/full.yaml
PYTHONPATH=src python -m neurothermo_per_cell.cli run --config configs/smoke.yaml
PYTHONPATH=src python -m neurothermo_per_cell.cli run --config configs/preliminary.yaml
PYTHONPATH=src python -m neurothermo_per_cell.cli run --config configs/full.yaml
```

Resume after interruption:

```bash
PYTHONPATH=src python -m neurothermo_per_cell.cli run \
  --config configs/preliminary.yaml --resume
```

Profile likelihood is deliberately not part of the preliminary run. After
reviewing `runs/preliminary/summary/full_cell_fits.csv`, place representative
converged cell IDs in `configs/profile.yaml` and run:

```bash
PYTHONPATH=src python -m neurothermo_per_cell.cli profile --config configs/profile.yaml
```

## Run sizes

| profile | cells | seeds | fits per cell/model | tasks |
|---|---:|---:|---:|---:|
| smoke | 2 | 1 | full + 1 CV | 8 |
| preliminary | 24 | 1 | full + 3 CV | 192 |
| full | 24 | 5 | full + 3 CV | 960 |

Preliminary and full profiles use four independent worker processes. Reduce
`workers` in the YAML file if the server has fewer available CPU cores. Do not
start two pipeline processes in the same output directory.

The full run must not be started merely because the preliminary run finishes.
Required stop/go evidence is: a high fraction of `CONVERGED` full-cell fits,
stable full-fit parameters across seeds in a targeted repeat, acceptable held-out
F–J/recruitment/latency errors, and informative profile likelihood for parameters
that will be interpreted biologically.

`parameter_stability_across_seeds.csv` is informative only when a profile uses
more than one seed. A single-seed preliminary run cannot establish parameter
stability; it only screens structural fit quality and identifies candidate cells
for a smaller multi-seed/profile calculation.

## Important limits

- Per-cell fitting removes between-cell averaging; it does not create information
  absent from F–J and first-spike-latency endpoints.
- Five dynamical parameters from 13 current levels can remain practically
  non-identifiable.
- A population claim requires comparing distributions of cell-level estimates,
  accounting for animals and parameter uncertainty. The pipeline does not treat
  13 sweeps as 13 independent cells.
- The equations and bounds are frozen in `models.py`; changing them requires a new
  version and output directory.
