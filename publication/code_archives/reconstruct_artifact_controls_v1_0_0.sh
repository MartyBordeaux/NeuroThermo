#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

cat artifact_controls_v1_0_0_source.tar.xz.b64.part* | tr -d '\n' | base64 -d > artifact_controls_v1_0_0_source.tar.xz

printf '%s  %s\n' \
  '8d1e6c99a4d25f89acd0e11c613830d4b6d86473835fd5346f8215c09ba3f988' \
  'artifact_controls_v1_0_0_source.tar.xz' | sha256sum -c -

mkdir -p ../code
tar -xJf artifact_controls_v1_0_0_source.tar.xz -C ../code

echo "Artifact-controls source extracted to publication/code/neurothermo_artifact_controls_v1_0_0/"
