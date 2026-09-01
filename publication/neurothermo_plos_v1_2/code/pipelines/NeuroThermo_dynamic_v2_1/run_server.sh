#!/bin/sh
set -eu
CONFIG="${1:-configs/server_dynamic_v2_1.yaml}"
python3 -m dynamic_v2.cli run --config "$CONFIG"
