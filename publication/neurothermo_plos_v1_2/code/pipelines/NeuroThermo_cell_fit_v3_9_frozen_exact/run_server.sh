#!/usr/bin/env bash
set -euo pipefail
CONFIG="${1:-configs/server_cellfit_v3_9.yaml}"
python3 -m hr_cell_fit.cli validate --config "$CONFIG"
python3 -m hr_cell_fit.cli run --config "$CONFIG"
echo "Fit complete. Review results_cellfit_v3_9/joint_fit_visual_audit_v3_9.pdf and edit review/comment in cell_fit_summary.csv before identify."
