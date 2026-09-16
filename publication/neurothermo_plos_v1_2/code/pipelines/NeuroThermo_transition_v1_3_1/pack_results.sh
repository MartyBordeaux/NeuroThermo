#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
DIR="${1:-results_transition_ensemble_v1_3}"
OUT="${2:-results_transition_ensemble_v1_3.zip}"
rm -f "$OUT"
zip -qr "$OUT" "$DIR"
echo "$HERE/$OUT"
