#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PACKAGE_DIR"
nohup ./run_server.sh "${1:-}" > kl_convergence_v1_0_1.log 2>&1 &
echo $! > kl_convergence_v1_0_1.pid
echo "PID $(cat kl_convergence_v1_0_1.pid)"
echo "Log: kl_convergence_v1_0_1.log"
