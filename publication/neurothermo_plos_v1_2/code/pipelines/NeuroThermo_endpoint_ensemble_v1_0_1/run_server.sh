#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
CONFIG="${1:-configs/server_endpoint_v1_0.yaml}"
python3 -m endpoint_v1.cli validate --config "$CONFIG"
python3 -m endpoint_v1.cli run --config "$CONFIG"
