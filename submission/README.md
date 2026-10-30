# Homework 2 submission

`report.pdf` is the final write-up. It contains the results, references,
Agent Audits, and the completed independent student checks.

## Reproduce the report

From the repository root, run:

```bash
bash submission/run_all.sh
```

This command creates the Python 3.12 environment, downloads and SHA-256
checks the inputs, runs every analysis in the listed order, and rebuilds
`submission/report.pdf`. It takes about four minutes on CPU. Tectonic is
downloaded to `.venv/bin` only if it is not already available.

## Analysis order

| Task | Scripts | Checkpoints |
|---|---|---|
| Single-cell PBMC3k | `single_cell/code/01_qc_baseline.py` through `04_ambiguous_population.py` | S1--S5 |
| Spatial mouse brain | `spatial/code/01_expression_domains.py` through `04_depth_confound_check.py` | T1--T5 |
| EGFR protein | `protein/code/01_map_and_measure.py` through `03_render_structure.py` | P1--P4 |
| Report | `report/build_report.sh` | PDF |

The numbered scripts run in the order above. Shared settings and seeds are in
each task's `code/common.py` (or `protein/code/structlib.py`): single-cell
seed 0 and subsample seed 42; spatial seed 0, with seed 1 for the secondary
null; spatial nulls use 1,000 permutations.

## Inputs

- PBMC3k matrix: course repository (`datasets/pbmc3k/`; 10x Genomics).
- Visium adult mouse brain matrix and spatial archive: course repository
  (`datasets/visium_mouse_brain/`; 10x Genomics).
- EGFR structures: course repository for 4HJO; RCSB for 1M17, 2ITX, 2ITN,
  2ITV, 2JIT, and 4I22.
- Protein annotations: UniProt P00533 and ERBB paralogs; AlphaFold DB
  AF-P00533-F1 v6, pLDDT, and AlphaMissense.

`fetch_data.sh` records and verifies SHA-256 checksums for the course inputs.

## What to submit

Submit this `submission/` directory, including `report.pdf`, this README,
the numbered scripts, and the tables and figures under each task's `outputs/`
folder. The `*.h5ad` files are regenerated intermediates (about 150 MB) and
may be left out of a zip file.
