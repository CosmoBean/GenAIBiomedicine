#!/usr/bin/env bash
# Reproduce every result in report.pdf, in order.
# Usage (from the repository root that contains submission/ and data/):
#   bash submission/run_all.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# 0. environment (Python 3.12; exact versions in submission/requirements.txt)
if [ ! -x .venv/bin/python ]; then
  uv venv --python 3.12 .venv
  uv pip install --python .venv/bin/python -r submission/requirements.txt
fi
PY="$ROOT/.venv/bin/python"

# 1. data (course repository copies; checksums verified against datasets/README.md)
bash submission/fetch_data.sh

# 2. Task 2 - single cell
cd "$ROOT/submission/single_cell/code"
$PY 01_qc_baseline.py
$PY 02_resolution_sweep.py
$PY 03_downsample.py
$PY 04_ambiguous_population.py

# 3. Task 3 - spatial
cd "$ROOT/submission/spatial/code"
$PY 01_expression_domains.py
$PY 02_spatial_statistics_and_null.py
$PY 03_svg_vs_de.py
$PY 04_depth_confound_check.py

# 4. Task 1 - protein
cd "$ROOT/submission/protein/code"
$PY 01_map_and_measure.py
$PY 02_mechanism_checks.py
$PY 03_render_structure.py

# 5. report
cd "$ROOT/submission/report"
$PY build_report.py
