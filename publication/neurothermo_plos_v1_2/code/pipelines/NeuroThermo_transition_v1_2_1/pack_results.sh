#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.."
DIR="results_transition_ensemble_v1_2_1"
ZIP="results_transition_ensemble_v1_2_1.zip"
rm -f "$ZIP"
zip -qr "$ZIP" "$DIR"
echo "$PWD/$ZIP"
