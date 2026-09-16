#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
RELEASE_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"
FROZEN_INPUT="$(cd "${1:-$RELEASE_ROOT/data/inputs}" && pwd)"
cd "$PACKAGE_DIR"
export PYTHONPATH="$PACKAGE_DIR${PYTHONPATH:+:$PYTHONPATH}"
export MPLCONFIGDIR="$PACKAGE_DIR/.matplotlib"
mkdir -p "$MPLCONFIGDIR"

python -m nonequilibrium_geometry validate \
  --config configs/server_nonequilibrium_geometry_v1_0_1.yaml \
  --frozen-dir "$FROZEN_INPUT"

python -m nonequilibrium_geometry run \
  --config configs/server_nonequilibrium_geometry_v1_0_1.yaml \
  --frozen-dir "$FROZEN_INPUT"
