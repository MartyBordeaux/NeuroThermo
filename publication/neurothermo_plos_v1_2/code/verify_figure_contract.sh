#!/usr/bin/env bash
set -euo pipefail

RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$RELEASE_ROOT/results/figures"
MODE="${1:-check}"

stems=(
  Fig1_endpoint_phenotype
  Fig2_transition_staging
  Fig3_intrinsic_drive_decomposition
  Fig4_thermodynamic_information_geometry
  Fig5_nonequilibrium_geometry
  FigS1_support_restricted_dynamics
  FigS2_multiseed_marker_robustness
  FigS3_KL_full_coverage_convergence
)

case "$MODE" in
  clean)
    find "$DIR" -maxdepth 1 -type f \( -name 'Fig*.pdf' -o -name 'Fig*.png' \) -delete
    ;;
  check)
    expected=$(mktemp)
    actual=$(mktemp)
    trap 'rm -f "$expected" "$actual"' EXIT
    for stem in "${stems[@]}"; do
      for ext in pdf png; do
        f="$DIR/$stem.$ext"
        test -s "$f"
        printf '%s\n' "$stem.$ext" >> "$expected"
      done
    done
    find "$DIR" -maxdepth 1 -type f \( -name 'Fig*.pdf' -o -name 'Fig*.png' \) -printf '%f\n' | sort > "$actual"
    sort -o "$expected" "$expected"
    diff -u "$expected" "$actual"
    test "$(wc -l < "$actual")" -eq 16
    echo FIGURE_CONTRACT_16_FILES_PASS
    ;;
  *)
    echo "Usage: $0 {clean|check}" >&2
    exit 2
    ;;
esac
