#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PIPE="$RELEASE_ROOT/code/pipelines/NeuroThermo_stage1_qc_fixed"
WT_ROOT="$RELEASE_ROOT/data/raw/WT"
SCA3_ROOT="$RELEASE_ROOT/data/raw/SCA3"
OUTPUT="${OUTPUT:-$RELEASE_ROOT/results/recomputed/stage1_qc_fixed}"

mkdir -p "$(dirname "$OUTPUT")"
cd "$PIPE"
python3 spike_qc_calibrated.py \
  --wt-root "$WT_ROOT" \
  --sca3-root "$SCA3_ROOT" \
  --config config.json \
  --qc-rules qc_rules.csv \
  --fixed-qc qc2.csv \
  --output "$OUTPUT"
