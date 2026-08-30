#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

case "${1:-}" in
  kl)
    exec "$RELEASE_ROOT/code/pipelines/NeuroThermo_KL_convergence_v1_0_1/run_server.sh"
    ;;
  nonequilibrium)
    exec "$RELEASE_ROOT/code/pipelines/NeuroThermo_nonequilibrium_geometry_v1_0_1/run_server.sh"
    ;;
  *)
    echo "Usage: $0 {kl|nonequilibrium}" >&2
    exit 2
    ;;
esac
