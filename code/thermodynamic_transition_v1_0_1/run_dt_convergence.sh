#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python3 -m thermo_v1_0_1.cli run --config configs/convergence_dt025_thermodynamic_v1_0_1.yaml
python3 -m thermo_v1_0_1.sensitivity dt \
  --reference results_thermodynamic_transition_v1_0_1 \
  --halfdt results_thermodynamic_transition_v1_0_1_dt025 \
  --out results_thermodynamic_transition_v1_0_1_dt_convergence
