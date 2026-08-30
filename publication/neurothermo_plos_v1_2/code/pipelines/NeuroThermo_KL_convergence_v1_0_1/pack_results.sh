#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR="results_kl_convergence_v1_0_1"
test -f "$RESULTS_DIR/KL_CONVERGENCE_VERDICT.json"
test -f "$RESULTS_DIR/RUN_SUMMARY.json"
zip -qr -FS neurothermo_kl_convergence_results_v1_0_1.zip "$RESULTS_DIR" kl_convergence_v1_0_1.log
unzip -t neurothermo_kl_convergence_results_v1_0_1.zip | tail -1
ls -lh neurothermo_kl_convergence_results_v1_0_1.zip
