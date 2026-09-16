#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PACKAGE_DIR"
export PYTHONPATH="$PACKAGE_DIR${PYTHONPATH:+:$PYTHONPATH}"
export MPLCONFIGDIR="$PACKAGE_DIR/.matplotlib"
mkdir -p "$MPLCONFIGDIR"

python -m unittest discover -s tests -v
python -m nonequilibrium_geometry validate \
  --config configs/smoke_nonequilibrium_geometry_v1_0_1.yaml \
  --frozen-dir frozen_smoke
python -m nonequilibrium_geometry run \
  --config configs/smoke_nonequilibrium_geometry_v1_0_1.yaml \
  --frozen-dir frozen_smoke

test -s results_smoke_nonequilibrium_geometry_v1_0_1/FORMALISM_VERDICT.json
test -s results_smoke_nonequilibrium_geometry_v1_0_1/ensemble_geometry.csv
test -s results_smoke_nonequilibrium_geometry_v1_0_1/preflight_endpoint_membership.csv
test -s results_smoke_nonequilibrium_geometry_v1_0_1/markov_cache/scenario_00000_seed_20260818.npz
test -s results_smoke_nonequilibrium_geometry_v1_0_1/figures/Fig_nonequilibrium_geometry.png
python - <<'PY'
import json
from pathlib import Path
import numpy as np
import pandas as pd

root = Path("results_smoke_nonequilibrium_geometry_v1_0_1")
geometry = pd.read_csv(root / "stationary_geometry.csv")
if not np.allclose(geometry["path_fi_xyz"], geometry["hs_force_variance"]):
    raise SystemExit("Centered FI/force-variance invariant failed")
fluctuation = pd.read_csv(root / "fluctuation_relations.csv")
if not (fluctuation["n_unique_path_positions"] == 5).all():
    raise SystemExit("Smoke protocol uniqueness failed")
cycles = pd.read_csv(root / "markov_cycle_affinities.csv")
required = {"scenario_id", "biological_pair_key", "seed", "p"}
if not required.issubset(cycles.columns):
    raise SystemExit("Cycle provenance failed")
verdict = json.loads((root / "FORMALISM_VERDICT.json").read_text())
if verdict["formalism_decision_basis"] != "coarse_markov_time_reversal":
    raise SystemExit("Formalism basis failed")
PY
echo "Smoke test PASS"
