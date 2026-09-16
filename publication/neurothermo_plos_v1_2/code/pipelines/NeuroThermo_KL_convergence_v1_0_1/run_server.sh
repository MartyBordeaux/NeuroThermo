#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
RELEASE_ROOT="$(cd "$PACKAGE_DIR/../../.." && pwd)"
FROZEN_DIR="${1:-$RELEASE_ROOT/data/inputs}"
cd "$PACKAGE_DIR"
export PYTHONPATH="$PACKAGE_DIR${PYTHONPATH:+:$PYTHONPATH}"
python3 -m kl_convergence validate --config configs/server_kl_convergence_v1_0_1.yaml --frozen-dir "$FROZEN_DIR"
python3 -m kl_convergence run --config configs/server_kl_convergence_v1_0_1.yaml --frozen-dir "$FROZEN_DIR"
