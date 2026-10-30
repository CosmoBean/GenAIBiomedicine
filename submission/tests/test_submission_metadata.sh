#!/usr/bin/env bash
# Fast consistency check for the grader-facing submission entry points.
set -euo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
report="$root/submission/report/report.tex"
readme="$root/submission/README.md"

rg -q 'My interpretation and validation' "$report"
! rg -q 'Student verification \(to complete before submission\)|needs student verification' "$report" "$readme"
rg -q 'report.pdf' "$readme"
rg -q 'bash submission/run_all.sh' "$readme"
