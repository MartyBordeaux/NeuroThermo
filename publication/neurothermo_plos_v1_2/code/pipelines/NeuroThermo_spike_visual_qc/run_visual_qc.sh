#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RELEASE_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

WT_ROOT="${WT_ROOT:-$RELEASE_ROOT/data/raw/WT}"
SCA3_ROOT="${SCA3_ROOT:-$RELEASE_ROOT/data/raw/SCA3}"
OUTPUT="${OUTPUT:-$RELEASE_ROOT/results/recomputed/spike_visual_qc}"
CONFIG="${CONFIG:-$SCRIPT_DIR/config.json}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/neurothermo_spike_visual_qc_mpl}"
mkdir -p "$MPLCONFIGDIR"

cd "$SCRIPT_DIR"
python3 spike_visual_qc.py \
  --wt-root "$WT_ROOT" \
  --sca3-root "$SCA3_ROOT" \
  --config "$CONFIG" \
  --output "$OUTPUT"
