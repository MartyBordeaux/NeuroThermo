#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CELLFIT="$RELEASE_ROOT/code/pipelines/NeuroThermo_cell_fit_v3_9_frozen_exact"
CELLFIT_CONFIG="$CELLFIT/configs/publication_cellfit_v3_9.yaml"

prepare() {
  python3 "$RELEASE_ROOT/code/prepare_calibration.py"
}

preflight() {
  prepare
  python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
}

usage() {
  cat >&2 <<'EOF'
Usage: code/run_full_analyses.sh {prepare|preflight|cellfit-validate|kl|nonequilibrium}

  prepare           Reconstruct and SHA-256 verify frozen calibration CSVs.
  preflight         Prepare inputs and run strict clean-clone release preflight.
  cellfit-validate  Validate exact frozen v3.9 cohort/input bundle; does not optimize.
  kl                Run the final KL convergence pipeline.
  nonequilibrium    Run the final nonequilibrium-geometry pipeline.
EOF
}

case "${1:-}" in
  prepare)
    prepare
    ;;
  preflight)
    preflight
    ;;
  cellfit-validate)
    preflight
    cd "$CELLFIT"
    python3 -m hr_cell_fit.cli validate --config "$CELLFIT_CONFIG"
    ;;
  kl)
    exec "$RELEASE_ROOT/code/pipelines/NeuroThermo_KL_convergence_v1_0_1/run_server.sh"
    ;;
  nonequilibrium)
    exec "$RELEASE_ROOT/code/pipelines/NeuroThermo_nonequilibrium_geometry_v1_0_1/run_server.sh"
    ;;
  *)
    usage
    exit 2
    ;;
esac
