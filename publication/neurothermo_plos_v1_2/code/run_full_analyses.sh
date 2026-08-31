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
TRANSITION121="$RELEASE_ROOT/code/pipelines/NeuroThermo_transition_v1_2_1"
TRANSITION121_CONFIG="$TRANSITION121/configs/publication_transition_v1_2_1.yaml"
TRANSITION121_FROZEN_CONFIG="$TRANSITION121/configs/publication_transition_v1_2_1_validate_frozen.yaml"
TRANSITION13="$RELEASE_ROOT/code/pipelines/NeuroThermo_transition_v1_3_frozen_exact"
TRANSITION13_CONFIG="$TRANSITION13/configs/publication_transition_v1_3.yaml"
TRANSITION13_FROZEN_CONFIG="$TRANSITION13/configs/publication_transition_v1_3_validate_frozen.yaml"

prepare_calibration() { python3 "$RELEASE_ROOT/code/prepare_calibration.py"; }
prepare_upstream() { python3 "$RELEASE_ROOT/code/prepare_upstream_inputs.py"; }
prepare_characterization() { prepare_upstream; python3 "$RELEASE_ROOT/code/prepare_characterization_input.py"; }
prepare_transition12() { python3 "$RELEASE_ROOT/code/prepare_transition_v1_2_inputs.py"; }

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
Usage: code/run_full_analyses.sh COMMAND

Core commands:
  prepare | prepare-upstream | preflight | raw-integrity | qc-tests | qc-recompute
  cellfit-validate | characterization | dynamic-validate | dynamic | endpoint-validate | endpoint
  transition-frozen | transition-v1-0-validate | transition-v1-0
  transition-v1-1-validate | transition-v1-1 | prepare-transition-v1-2
  transition-v1-2-validate | transition-v1-2
  transition-v1-2-1-validate-frozen | transition-v1-2-1-validate | transition-v1-2-1
  transition-v1-3-validate-frozen | transition-v1-3-validate | transition-v1-3
  transition-integrity | figure-source | figures-python | kl | nonequilibrium
EOF
}

case "${1:-}" in
  prepare) prepare_calibration ;;
  prepare-upstream) prepare_upstream ;;
  preflight) preflight ;;
  raw-integrity) python3 "$RELEASE_ROOT/code/validate_raw_data.py" ;;
  qc-tests) cd "$QC"; python3 -m pytest -q ;;
  qc-recompute) exec bash "$RELEASE_ROOT/code/run_qc_publication.sh" ;;
  cellfit-validate)
    prepare_calibration
    python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
    cd "$CELLFIT"; python3 -m hr_cell_fit.cli validate --config "$CELLFIT_CONFIG"
    ;;
  characterization)
    prepare_characterization
    python3 "$CHAR/run_characterization.py" \
      --results "$RELEASE_ROOT/data/v3_9_results_full" \
      --animal-map "$RELEASE_ROOT/results/recomputed/characterization_inputs/NeuroThermo_animal_id_recovery.xlsx" \
      --outdir "$RELEASE_ROOT/results/recomputed/characterization_v1_0"
    ;;
  dynamic-validate) prepare_upstream; cd "$DYNAMIC"; python3 -m dynamic_v2.cli validate --config "$DYNAMIC_CONFIG" ;;
  dynamic) prepare_upstream; cd "$DYNAMIC"; python3 -m dynamic_v2.cli run --config "$DYNAMIC_CONFIG" ;;
  endpoint-validate) cd "$ENDPOINT"; python3 -m endpoint_v1.cli validate --config "$ENDPOINT_CONFIG" ;;
  endpoint) cd "$ENDPOINT"; python3 -m endpoint_v1.cli run --config "$ENDPOINT_CONFIG" ;;
  transition-frozen) python3 "$RELEASE_ROOT/code/validate_transition_frozen.py" ;;
  transition-v1-0-validate) cd "$TRANSITION10"; python3 -m transition_v1.cli validate --config "$TRANSITION10_CONFIG" ;;
  transition-v1-0) cd "$TRANSITION10"; python3 -m transition_v1.cli run --config "$TRANSITION10_CONFIG" ;;
  transition-v1-1-validate) cd "$TRANSITION11"; python3 -m transition_v1_1.cli validate --config "$TRANSITION11_CONFIG" ;;
  transition-v1-1) cd "$TRANSITION11"; python3 -m transition_v1_1.cli run --config "$TRANSITION11_CONFIG" ;;
  prepare-transition-v1-2) prepare_transition12 ;;
  transition-v1-2-validate) prepare_transition12; cd "$TRANSITION12"; python3 -m transition_v1_2.cli validate --config "$TRANSITION12_CONFIG" ;;
  transition-v1-2) prepare_transition12; cd "$TRANSITION12"; python3 -m transition_v1_2.cli run --config "$TRANSITION12_CONFIG" ;;
  transition-v1-2-1-validate-frozen) cd "$TRANSITION121"; python3 -m transition_v1_2_1.cli validate --config "$TRANSITION121_FROZEN_CONFIG" ;;
  transition-v1-2-1-validate) cd "$TRANSITION121"; python3 -m transition_v1_2_1.cli validate --config "$TRANSITION121_CONFIG" ;;
  transition-v1-2-1) cd "$TRANSITION121"; python3 -m transition_v1_2_1.cli run --config "$TRANSITION121_CONFIG" ;;
  transition-v1-3-validate-frozen) cd "$TRANSITION13"; python3 -m transition_v1_3.cli validate --config "$TRANSITION13_FROZEN_CONFIG" ;;
  transition-v1-3-validate) cd "$TRANSITION13"; python3 -m transition_v1_3.cli validate --config "$TRANSITION13_CONFIG" ;;
  transition-v1-3) cd "$TRANSITION13"; python3 -m transition_v1_3.cli run --config "$TRANSITION13_CONFIG" ;;
  transition-integrity) python3 "$RELEASE_ROOT/code/validate_transition_results.py" ;;
  figure-source) python3 "$RELEASE_ROOT/code/assemble_figure_source.py" ;;
  figures-python) python3 "$RELEASE_ROOT/code/figures/python/render_figures.py"; python3 "$RELEASE_ROOT/code/figures/python/render_nonequilibrium_summary.py" ;;
  kl) exec "$RELEASE_ROOT/code/pipelines/NeuroThermo_KL_convergence_v1_0_1/run_server.sh" ;;
  nonequilibrium) exec "$RELEASE_ROOT/code/pipelines/NeuroThermo_nonequilibrium_geometry_v1_0_1/run_server.sh" ;;
  *) usage; exit 2 ;;
esac
