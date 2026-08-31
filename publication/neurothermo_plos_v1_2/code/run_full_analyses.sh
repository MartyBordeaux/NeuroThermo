#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

usage() {
  cat >&2 <<'EOF'
Usage: code/run_full_analyses.sh {preflight|upstream-status|kl|nonequilibrium}

  preflight        Run strict clean-clone release preflight.
  upstream-status  Show upstream readiness without forcing a non-zero exit.
  kl               Run the final KL convergence pipeline.
  nonequilibrium   Run the final nonequilibrium-geometry pipeline.

The upstream cell-fit -> transition chain is intentionally not launched until
all provenance-critical calibration/frozen inputs pass preflight. This avoids
silent substitution of historical server files or machine-specific paths.
EOF
}

case "${1:-}" in
  preflight)
    exec python3 "$RELEASE_ROOT/code/preflight_release.py" --strict
    ;;
  upstream-status)
    exec python3 "$RELEASE_ROOT/code/preflight_release.py"
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
