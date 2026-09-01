#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RELEASE_ROOT"

# Rebuild every publication source table derived deterministically from frozen
# upstream layers and require numerical equivalence to immutable references.
python code/assemble_figure_source.py

# Render the complete frozen publication figure contract from one canonical
# entrypoint: 8 figure stems, each as PDF and PNG.
bash code/run_full_analyses.sh figures-python

echo "Publication source-data lineage verified and all canonical figures written to $RELEASE_ROOT/results/figures"
