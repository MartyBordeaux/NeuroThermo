#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
cat publication_source_code.tar.xz.b64.part* | tr -d '\n' | base64 -d > publication_source_code.tar.xz
sha256sum -c publication_source_code.sha256
mkdir -p ../code
tar -xJf publication_source_code.tar.xz -C ../code
echo "Source code extracted to publication/code/"
