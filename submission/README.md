# Homework 2: submission guide

## Start here

Open **`report.pdf`** first. It is the complete report: results, figures,
references, Agent Audits, and my interpretation and validation checks.

## Run the analysis (optional)

From the repository root, run:

```bash
bash submission/run_all.sh
```

This creates the Python environment, downloads the required inputs, runs all
analyses, and rebuilds `submission/report.pdf`. It takes about four minutes on
a CPU. Python 3.12 and the exact package versions are listed in
`requirements.txt`.

## Where to find things

| If you want to see... | Open... |
|---|---|
| The final answers | `report.pdf` |
| Protein analysis (P1--P4) | `protein/code/` and `protein/outputs/` |
| Single-cell analysis (S1--S5) | `single_cell/code/` and `single_cell/outputs/` |
| Spatial analysis (T1--T5) | `spatial/code/` and `spatial/outputs/` |
| The report source | `report/` |

The numbered scripts run in order. Settings and random seeds are recorded in
each task's `code/common.py` file (or `protein/code/structlib.py`).

## Contents of this folder

- `report.pdf` — final report.
- `README.md`, `requirements.txt`, `fetch_data.sh`, and `run_all.sh` — setup
  and reproduction instructions.
- `protein/`, `single_cell/`, and `spatial/` — code plus the figures and
  tables used in the report.

The large `*.h5ad` files are temporary, regenerated analysis objects. They
are intentionally not included in the submission ZIP.
