#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PACKAGE_DIR"
nohup ./run_server.sh "${1:-}" > nonequilibrium_geometry_v1_0_1.log 2>&1 &
PID=$!
echo "$PID" > nonequilibrium_geometry_v1_0_1.pid
echo "Started PID $PID"
echo "Log: $PACKAGE_DIR/nonequilibrium_geometry_v1_0_1.log"
