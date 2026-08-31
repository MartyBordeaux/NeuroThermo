#!/usr/bin/env bash
set -euo pipefail

WT_ROOT="${WT_ROOT:-/root/neurothermo/WT}"
SCA3_ROOT="${SCA3_ROOT:-/root/neurothermo/SCA3}"
OUTPUT="${OUTPUT:-results_spike_visual_qc}"
CONFIG="${CONFIG:-config.json}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/neurothermo_spike_visual_qc_mpl}"
mkdir -p "$MPLCONFIGDIR"

python3 spike_visual_qc.py \
  --wt-root "$WT_ROOT" \
  --sca3-root "$SCA3_ROOT" \
  --config "$CONFIG" \
  --output "$OUTPUT"
