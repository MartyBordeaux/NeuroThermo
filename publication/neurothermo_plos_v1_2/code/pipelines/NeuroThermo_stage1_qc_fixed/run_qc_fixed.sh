#!/usr/bin/env bash
set -euo pipefail

WT_ROOT="${WT_ROOT:-/root/neurothermo/WT}"
SCA3_ROOT="${SCA3_ROOT:-/root/neurothermo/SCA3}"
OUTPUT="${OUTPUT:-results_stage1_qc_fixed}"

python3 spike_qc_calibrated.py \
  --wt-root "$WT_ROOT" \
  --sca3-root "$SCA3_ROOT" \
  --config config.json \
  --qc-rules qc_rules.csv \
  --fixed-qc qc2.csv \
  --output "$OUTPUT"
