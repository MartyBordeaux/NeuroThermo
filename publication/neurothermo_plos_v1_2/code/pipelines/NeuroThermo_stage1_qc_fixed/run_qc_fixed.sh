#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RELEASE_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

WT_ROOT="${WT_ROOT:-$RELEASE_ROOT/data/raw/WT}"
SCA3_ROOT="${SCA3_ROOT:-$RELEASE_ROOT/data/raw/SCA3}"
OUTPUT="${OUTPUT:-$RELEASE_ROOT/results/recomputed/stage1_qc_fixed}"

cd "$SCRIPT_DIR"
python3 spike_qc_calibrated.py \
  --wt-root "$WT_ROOT" \
  --sca3-root "$SCA3_ROOT" \
  --config config.json \
  --qc-rules qc_rules.csv \
  --fixed-qc qc2.csv \
  --output "$OUTPUT"
