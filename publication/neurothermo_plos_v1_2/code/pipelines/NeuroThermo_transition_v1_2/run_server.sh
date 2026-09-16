#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
CFG="${1:-configs/server_transition_v1_2.yaml}"
python3 -m transition_v1_2.cli validate --config "$CFG"
python3 -m transition_v1_2.cli run --config "$CFG"
