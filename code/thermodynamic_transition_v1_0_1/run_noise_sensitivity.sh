#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python3 -m thermo_v1_0_1.cli run --config configs/server_thermodynamic_v1_0_1_noise_half.yaml
python3 -m thermo_v1_0_1.cli run --config configs/server_thermodynamic_v1_0_1_noise_double.yaml
python3 -m thermo_v1_0_1.sensitivity noise \
  --reference results_thermodynamic_transition_v1_0_1 \
  --half results_thermodynamic_transition_v1_0_1_noise_half \
  --double results_thermodynamic_transition_v1_0_1_noise_double \
  --out results_thermodynamic_transition_v1_0_1_noise_sensitivity
