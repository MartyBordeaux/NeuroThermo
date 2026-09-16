#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="${1:-results_transition_ensemble_v1_2}"
if [ ! -d "$OUT" ]; then
  echo "Result directory not found: $OUT" >&2
  exit 1
fi
rm -f "${OUT}.zip"
zip -qr "${OUT}.zip" "$OUT"
echo "Created: ${OUT}.zip"
