#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
CONFIG="${1:-configs/server_transition_v1_1.yaml}"
python3 -m transition_v1_1.cli validate --config "$CONFIG"
python3 -m transition_v1_1.cli run --config "$CONFIG"
