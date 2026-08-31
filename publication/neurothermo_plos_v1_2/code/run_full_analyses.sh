#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CELLFIT="$RELEASE_ROOT/code/pipelines/NeuroThermo_cell_fit_v3_9_frozen_exact"
CELLFIT_CONFIG="$CELLFIT/configs/publication_cellfit_v3_9.yaml"
DYNAMIC="$RELEASE_ROOT/code/pipelines/NeuroThermo_dynamic_v2_1"
DYNAMIC_CONFIG="$DYNAMIC/configs/publication_dynamic_v2_1.yaml"
ENDPOINT="$RELEASE_ROOT/code/pipelines/NeuroThermo_endpoint_ensemble_v1_0_1"
ENDPOINT_CONFIG="$ENDPOINT/configs/publication_endpoint_v1_0.yaml"

prepare_calibration() {
  python3 "$RELEASE_ROOT/code/prepare_calibration.py"
}

prepare_upstream() {
  python3 "$RELEASE_ROOT/code/prepare_upstream_inputs.py"
}

preflight() {
  prepare_calibration
  prepare_upstream
  python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
}

usage() {
  cat >&2 <<'EOF'
Usage: code/run_full_analyses.sh {prepare|prepare-upstream|preflight|cellfit-validate|dynamic-validate|dynamic|endpoint-validate|endpoint|transition-integrity|kl|nonequilibrium}

  prepare              Reconstruct and SHA-256 verify frozen calibration CSVs.
  prepare-upstream     Extract and SHA-256 verify imported full v3.9 results and dynamic-v2.1 frozen inputs.
  preflight            Prepare all frozen inputs and run strict clean-clone release preflight.
  cellfit-validate     Validate exact frozen v3.9 cohort/input bundle; does not optimize.
  dynamic-validate     Validate experimental-support-restricted dynamic-v2.1 frozen input layer.
  dynamic              Recompute dynamic-v2.1 from the frozen publication input layer.
  endpoint-validate    Validate endpoint-v1.0.1 frozen input layer.
  endpoint             Recompute endpoint-v1.0.1 from the frozen publication input layer.
  transition-integrity Validate cross-stage v1.0/v1.1/v1.2/v1.2.1/v1.3 frozen results and embedded input hashes.
  kl                   Run the final KL convergence pipeline.
  nonequilibrium       Run the final nonequilibrium-geometry pipeline.
EOF
}

case "${1:-}" in
  prepare)
    prepare_calibration
    ;;
  prepare-upstream)
    prepare_upstream
    ;;
  preflight)
    preflight
    ;;
  cellfit-validate)
    prepare_calibration
    python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
    cd "$CELLFIT"
    python3 -m hr_cell_fit.cli validate --config "$CELLFIT_CONFIG"
    ;;
  dynamic-validate)
    prepare_upstream
    cd "$DYNAMIC"
    python3 -m dynamic_v2.cli validate --config "$DYNAMIC_CONFIG"
    ;;
  dynamic)
    prepare_upstream
    cd "$DYNAMIC"
    python3 -m dynamic_v2.cli run --config "$DYNAMIC_CONFIG"
    ;;
  endpoint-validate)
    cd "$ENDPOINT"
    python3 -m endpoint_v1.cli validate --config "$ENDPOINT_CONFIG"
    ;;
  endpoint)
    cd "$ENDPOINT"
    python3 -m endpoint_v1.cli run --config "$ENDPOINT_CONFIG"
    ;;
  transition-integrity)
    python3 "$RELEASE_ROOT/code/validate_transition_results.py"
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
