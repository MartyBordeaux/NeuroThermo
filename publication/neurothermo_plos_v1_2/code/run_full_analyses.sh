#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
QC="$RELEASE_ROOT/code/pipelines/NeuroThermo_stage1_qc_fixed"
CELLFIT="$RELEASE_ROOT/code/pipelines/NeuroThermo_cell_fit_v3_9_frozen_exact"
CELLFIT_CONFIG="$CELLFIT/configs/publication_cellfit_v3_9.yaml"
CHAR="$RELEASE_ROOT/code/pipelines/NeuroThermo_characterization_v1_0"
DYNAMIC="$RELEASE_ROOT/code/pipelines/NeuroThermo_dynamic_v2_1"
DYNAMIC_CONFIG="$DYNAMIC/configs/publication_dynamic_v2_1.yaml"
ENDPOINT="$RELEASE_ROOT/code/pipelines/NeuroThermo_endpoint_ensemble_v1_0_1"
ENDPOINT_CONFIG="$ENDPOINT/configs/publication_endpoint_v1_0.yaml"
TRANSITION10="$RELEASE_ROOT/code/pipelines/NeuroThermo_transition_v1_0"
TRANSITION10_CONFIG="$TRANSITION10/configs/publication_transition_v1_0.yaml"
TRANSITION11="$RELEASE_ROOT/code/pipelines/NeuroThermo_transition_v1_1"
TRANSITION11_CONFIG="$TRANSITION11/configs/publication_transition_v1_1.yaml"
TRANSITION12="$RELEASE_ROOT/code/pipelines/NeuroThermo_transition_v1_2"
TRANSITION12_CONFIG="$TRANSITION12/configs/publication_transition_v1_2.yaml"

prepare_calibration() {
  python3 "$RELEASE_ROOT/code/prepare_calibration.py"
}

prepare_upstream() {
  python3 "$RELEASE_ROOT/code/prepare_upstream_inputs.py"
}

prepare_characterization() {
  prepare_upstream
  python3 "$RELEASE_ROOT/code/prepare_characterization_input.py"
}

prepare_transition12() {
  python3 "$RELEASE_ROOT/code/prepare_transition_v1_2_inputs.py"
}

preflight() {
  prepare_calibration
  prepare_upstream
  python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
  python3 "$RELEASE_ROOT/code/validate_raw_data.py"
  python3 "$RELEASE_ROOT/code/validate_transition_frozen.py"
  python3 "$RELEASE_ROOT/code/validate_transition_results.py"
}

usage() {
  cat >&2 <<'EOF'
Usage: code/run_full_analyses.sh {prepare|prepare-upstream|preflight|raw-integrity|qc-tests|qc-recompute|cellfit-validate|characterization|dynamic-validate|dynamic|endpoint-validate|endpoint|transition-frozen|transition-v1-0-validate|transition-v1-0|transition-v1-1-validate|transition-v1-1|prepare-transition-v1-2|transition-v1-2-validate|transition-v1-2|transition-integrity|kl|nonequilibrium}

  prepare                  Reconstruct and SHA-256 verify frozen calibration CSVs.
  prepare-upstream         Extract and SHA-256 verify imported full v3.9 results and dynamic-v2.1 frozen inputs.
  preflight                Run all static, raw-data, frozen-input and transition-chain integrity checks.
  raw-integrity            Verify all 50 raw ABF recordings against RAW_DATA_MANIFEST.tsv.
  qc-tests                 Run tests for the final stage1 fixed-QC implementation.
  qc-recompute             Recompute raw ABF -> candidate events -> fixed manual QC using release-relative paths.
  cellfit-validate         Validate exact frozen v3.9 cohort/input bundle; does not optimize.
  characterization         Recompute post-fit characterization from full v3.9 results and canonical animal map.
  dynamic-validate         Validate experimental-support-restricted dynamic-v2.1 frozen input layer.
  dynamic                  Recompute dynamic-v2.1 from the frozen publication input layer.
  endpoint-validate        Validate endpoint-v1.0.1 frozen input layer.
  endpoint                 Recompute endpoint-v1.0.1 from the frozen publication input layer.
  transition-frozen        Verify exact v1.0/v1.1/v1.2/v1.3 transition frozen inputs and SHA-256.
  transition-v1-0-validate Validate portable transition-v1.0 config and all frozen inputs.
  transition-v1-0          Recompute full transition-v1.0 ensemble; checkpoint/resume enabled.
  transition-v1-1-validate Validate v1.1 reprojection against imported v1.0 results and exact q75 reference.
  transition-v1-1          Recompute zero-simulation transition-v1.1 reprojection.
  prepare-transition-v1-2 Assemble v1.2 inputs from recomputed v1.1 and verify against historical frozen references.
  transition-v1-2-validate Validate portable transition-v1.2 config using assembled inputs.
  transition-v1-2          Recompute full 31x31 transition-v1.2 surface; checkpoint/resume enabled.
  transition-integrity     Validate cross-stage v1.0/v1.1/v1.2/v1.2.1/v1.3 frozen results and embedded input hashes.
  kl                       Run the final KL convergence pipeline.
  nonequilibrium           Run the final nonequilibrium-geometry pipeline.
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
  raw-integrity)
    python3 "$RELEASE_ROOT/code/validate_raw_data.py"
    ;;
  qc-tests)
    cd "$QC"
    python3 -m pytest -q
    ;;
  qc-recompute)
    exec bash "$RELEASE_ROOT/code/run_qc_publication.sh"
    ;;
  cellfit-validate)
    prepare_calibration
    python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
    cd "$CELLFIT"
    python3 -m hr_cell_fit.cli validate --config "$CELLFIT_CONFIG"
    ;;
  characterization)
    prepare_characterization
    python3 "$CHAR/run_characterization.py" \
      --results "$RELEASE_ROOT/data/v3_9_results_full" \
      --animal-map "$RELEASE_ROOT/results/recomputed/characterization_inputs/NeuroThermo_animal_id_recovery.xlsx" \
      --outdir "$RELEASE_ROOT/results/recomputed/characterization_v1_0"
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
  transition-frozen)
    python3 "$RELEASE_ROOT/code/validate_transition_frozen.py"
    ;;
  transition-v1-0-validate)
    cd "$TRANSITION10"
    python3 -m transition_v1.cli validate --config "$TRANSITION10_CONFIG"
    ;;
  transition-v1-0)
    cd "$TRANSITION10"
    python3 -m transition_v1.cli run --config "$TRANSITION10_CONFIG"
    ;;
  transition-v1-1-validate)
    cd "$TRANSITION11"
    python3 -m transition_v1_1.cli validate --config "$TRANSITION11_CONFIG"
    ;;
  transition-v1-1)
    cd "$TRANSITION11"
    python3 -m transition_v1_1.cli run --config "$TRANSITION11_CONFIG"
    ;;
  prepare-transition-v1-2)
    prepare_transition12
    ;;
  transition-v1-2-validate)
    prepare_transition12
    cd "$TRANSITION12"
    python3 -m transition_v1_2.cli validate --config "$TRANSITION12_CONFIG"
    ;;
  transition-v1-2)
    prepare_transition12
    cd "$TRANSITION12"
    python3 -m transition_v1_2.cli run --config "$TRANSITION12_CONFIG"
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
