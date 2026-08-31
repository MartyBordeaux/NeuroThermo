#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RELEASE_ROOT"

# Rebuild every Fig1--3 source table consumed by the canonical renderer and
# require numerical equivalence to the immutable publication references.
python code/assemble_figure_source.py

# Render every canonical main/supporting figure from committed publication
# source tables or frozen downstream result layers.
python code/figures/python/render_figures.py
python code/figures/python/render_fig4_multiseed.py
python code/figures/python/render_nonequilibrium_summary.py
python code/figures/python/render_supporting_robustness.py

echo "Publication source-data lineage verified and figures written to $RELEASE_ROOT/results/figures"
