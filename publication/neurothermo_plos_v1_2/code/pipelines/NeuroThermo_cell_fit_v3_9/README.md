# NeuroThermo HR cell fit v3.9

v3.9 is a one-factor global boundary-stress experiment built directly on the accepted v3.8 formulation. The only scientific change is extension of the lower bound of `s` from 0.25 to 0.05 for every cell.

Frozen formulation:

- 20 accepted cells;
- 18 primary multi-sweep cells and 2 single-sweep cells;
- 113 frozen spiking sweeps, 4884 selected spikes;
- fitted parameters: `b`, `r`, `s`, `kappa_I`;
- exact additive first-spike alignment; no time rescaling and no last-spike anchor;
- binary rheobase bracket;
- no explicit latency, ISI, adaptation, voltage or plateau loss.

`SCA3_05` remains an ordinary member of the accepted fit cohort. v3.9 contains no special rescue or sensitivity branch for that cell.

## Parameter space

```text
b       0.5       ... 7.0       linear
r       0.0001    ... 0.10      log
s       0.05      ... 15.0      linear   # only changed bound
kappa_I 0.0002    ... 2.0       log
```

All cells are refitted with the same strong global search. Final v3.8 solutions are supplied only as optimizer seeds and explicit baselines.

## Purpose

The experiment asks whether v3.8 estimates of `s`, especially WT values near 0.25, were artificially constrained by the lower bound. `s_boundary_stress_summary.csv` reports the v3.8 and v3.9 estimates side by side, whether the new optimum moved below 0.25, and its distance from the new lower bound.

Identifiability is evaluated in the v3.9 search domain, but the required separated-alternative distance is still defined by the original v3.6 reference ranges. Thus lowering `s_min` does not make an `IDENTIFIABLE` call artificially easier.

## Usage

```bash
cd ~/neurothermo/v4
unzip hr_cell_fit_v3.9.0.zip
cd hr_cell_fit_v3_9
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -e .

python3 -m hr_cell_fit.cli validate --config configs/server_cellfit_v3_9.yaml
./run_server.sh configs/server_cellfit_v3_9.yaml
```

Review:

- `results_cellfit_v3_9/joint_fit_visual_audit_v3_9.pdf`
- `results_cellfit_v3_9/cell_fit_summary.csv`
- `results_cellfit_v3_9/s_boundary_stress_summary.csv`
- `results_cellfit_v3_9/latency_alignment_summary.csv`
- `results_cellfit_v3_9/threshold_constraint_summary.csv`

After visual review:

```bash
python3 -m hr_cell_fit.cli identify --config configs/server_cellfit_v3_9.yaml
```
