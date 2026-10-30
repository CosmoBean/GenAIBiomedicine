#!/usr/bin/env bash
# Compile report.tex -> ../report.pdf with Tectonic (installed into .venv/bin on first use).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
TT="$HERE/../../.venv/bin/tectonic"
if [ ! -x "$TT" ]; then
  (cd "$(dirname "$TT")" && curl -sSL --proto '=https' https://drop-sh.fullyjustified.net | sh)
fi
cd "$HERE"
"$TT" report.tex
mv report.pdf ../report.pdf
echo "wrote $(cd .. && pwd)/report.pdf"
