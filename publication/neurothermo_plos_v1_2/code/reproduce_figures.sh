#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RELEASE_ROOT"

python code/figures/python/render_figures.py
python code/figures/python/render_nonequilibrium_summary.py
python code/figures/python/render_supporting_robustness.py

echo "Publication figures written to $RELEASE_ROOT/results/figures"
