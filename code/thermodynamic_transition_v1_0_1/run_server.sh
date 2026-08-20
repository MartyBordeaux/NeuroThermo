#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
CONFIG="${1:-configs/server_thermodynamic_v1_0_1.yaml}"
python3 -m thermo_v1_0_1.cli validate --config "$CONFIG"
python3 -m thermo_v1_0_1.cli run --config "$CONFIG"
