#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if command -v bibtex >/dev/null 2>&1; then BIB=bibtex
elif command -v bibtex.original >/dev/null 2>&1; then BIB=bibtex.original
else echo 'BibTeX is required.' >&2; exit 1; fi
for doc in main supplementary response_to_reviewers; do
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" > "$doc.build.log" 2>&1
  "$BIB" "$doc" >> "$doc.build.log" 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" >> "$doc.build.log" 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" >> "$doc.build.log" 2>&1
done
