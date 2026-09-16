#!/usr/bin/env bash
set -euo pipefail

python3 -m kl_convergence validate --config configs/smoke_kl_convergence_v1_0_1.yaml --frozen-dir frozen_smoke
python3 -m kl_convergence run --config configs/smoke_kl_convergence_v1_0_1.yaml --frozen-dir frozen_smoke
test -f smoke_results_kl_convergence_v1_0_1/KL_CONVERGENCE_VERDICT.json
test -f smoke_results_kl_convergence_v1_0_1/GRID_COVERAGE_AUDIT.json
test -f smoke_results_kl_convergence_v1_0_1/Fig_KL_convergence.pdf
python3 -c 'import json; p=json.load(open("smoke_results_kl_convergence_v1_0_1/GRID_COVERAGE_AUDIT.json")); assert p["all_tasks_full_coverage"] and p["minimum_retained_mass"] == 1.0'
python3 -m unittest discover -s tests -v
