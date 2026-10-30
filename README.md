# GenAI for Biomedicine — Homework 2

Scientific analysis with an AI agent workspace (CMU 11-781 / 02-741, Fall 2026):
EGFR variant structure analysis, PBMC3k single-cell robustness, and Visium spatial domains.

- Report: [`submission/report.pdf`](submission/report.pdf)
- Details, run order and data sources: [`submission/README.md`](submission/README.md)

Reproduce everything (CPU, ~4 min; needs `uv`):

```bash
bash submission/run_all.sh
```

This creates `.venv`, downloads and checksums the inputs into `data/`, runs all analyses and rebuilds the report.
