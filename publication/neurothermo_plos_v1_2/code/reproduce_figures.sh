#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RELEASE_ROOT"

# Rebuild deterministic publication source tables from frozen upstream layers.
python code/assemble_figure_source.py

# Remove every committed publication render so the replay cannot pass by
# inheriting an old PDF/PNG from the clone.
bash code/verify_figure_contract.sh clean

# Render the complete frozen publication figure contract from one canonical
# entrypoint: 8 figure stems, each as PDF and PNG.
bash code/run_full_analyses.sh figures-python

# Require exactly the 16 frozen release files and reject unexpected renders.
bash code/verify_figure_contract.sh check

echo "Publication source-data lineage verified and all 16 canonical figure files written to $RELEASE_ROOT/results/figures"
