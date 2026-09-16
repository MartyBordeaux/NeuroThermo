#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
CFG="${1:-configs/server_transition_v1_0.yaml}"
python3 -m transition_v1.cli validate --config "$CFG"
python3 -m transition_v1.cli run --config "$CFG"
