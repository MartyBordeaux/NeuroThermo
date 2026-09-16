#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

for package in \
  "$RELEASE_ROOT/code/pipelines/NeuroThermo_KL_convergence_v1_0_1" \
  "$RELEASE_ROOT/code/pipelines/NeuroThermo_nonequilibrium_geometry_v1_0_1"
do
  (cd "$package" && ./smoke_test.sh)
done
