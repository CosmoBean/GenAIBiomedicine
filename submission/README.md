# Homework 2: Scientific Analysis with AI Agent Workspaces

CMU 11-781 / 02-741, Fall 2026. The full write-up, with all three Agent Audits and the checkpoint evidence index, is in **`report.pdf`**.

The agent workspace was Claude Code with the Claude Opus 5.5 model. It ran on CPU only; no GPU or AlphaFold run was needed.

## Run everything

From the directory that contains `submission/`:

```bash
bash submission/run_all.sh     # about 4 min on CPU
```

`run_all.sh` does the following:
1. Creates `.venv` with uv from `requirements.txt` (Python 3.12).
2. Runs `fetch_data.sh` to download and SHA-256-check all inputs into `data/`.
3. Runs the scripts below in order.
4. Compiles `report.pdf` from LaTeX (Tectonic is downloaded into `.venv/bin` on first use).

| Order | Script | Checkpoints | Main outputs |
|---|---|---|---|
| 1 | `single_cell/code/01_qc_baseline.py` | S1, S2 | QC figures, baseline UMAP, Wilcoxon markers, annotation |
| 2 | `single_cell/code/02_resolution_sweep.py` | S3 | ARI/NMI, splits and merges, label stability, seed control |
| 3 | `single_cell/code/03_downsample.py` | S4 | 70 % subsample (seed 42) comparison, extra seeds and resolutions |
| 4 | `single_cell/code/04_ambiguous_population.py` | S5 | naive vs memory CD4 analysis |
| 5 | `spatial/code/01_expression_domains.py` | T1 | expression-only domains (coordinates asserted absent) |
| 6 | `spatial/code/02_spatial_statistics_and_null.py` | T2, T5 | hex-grid graph, neighbourhood enrichment, Moran's I, 1,000-repeat label-permutation null + graph-rebuild check |
| 7 | `spatial/code/03_svg_vs_de.py` | T3, T4 | Moran's I vs Wilcoxon DE (z and log2FC) on the same 200 HVGs; depth/mito covariates |
| 8 | `spatial/code/04_depth_confound_check.py` | T5 | re-clustering with depth and mito regressed out |
| 9 | `protein/code/01_map_and_measure.py` | P1, P2 | sequence check, SEQRES mapping, distances, SASA, pLDDT, paralogs |
| 10 | `protein/code/02_mechanism_checks.py` | P2, P3 | T790M superpositions, G719S P-loop, L858 environment |
| 11 | `protein/code/03_render_structure.py` | P4 | annotated 4HJO figure plus interactive HTML |
| 12 | `report/build_report.sh` | R | `report.pdf` from `report/report.tex` (LaTeX, compiled with Tectonic) |

Shared settings and seeds are in `single_cell/code/common.py`, `spatial/code/common.py` and `protein/code/structlib.py`.
- **Single cell:** seed 0; subsample seed 42.
- **Spatial:** seed 0 for clustering and Null A; seed 1 for Null B; 1,000 permutations.

## Data sources

| Task | File | Source |
|---|---|---|
| 2 | `pbmc3k_filtered_gene_bc_matrices.tar.gz` | course repo `datasets/pbmc3k/` (10x Genomics PBMC3k) |
| 3 | `V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5`, `V1_Adult_Mouse_Brain_spatial.tar.gz` | course repo `datasets/visium_mouse_brain/` (10x Genomics) |
| 1 | `4hjo.pdb` | course repo `datasets/protein/` (RCSB 4HJO) |
| 1 | 1M17, 2ITX, 2ITN, 2ITV, 2JIT, 4I22 | files.rcsb.org |
| 1 | UniProt P00533 (+ paralogs P04626, P21860, Q15303) | rest.uniprot.org |
| 1 | AF-P00533-F1 v6 model, pLDDT, AlphaMissense | alphafold.ebi.ac.uk |

## Folder layout

```
submission/
├── report.pdf            final report
├── report/               report.tex + cisXXX.cls (UPenn-style HW class) + build_report.sh
├── README.md  run_all.sh  fetch_data.sh  requirements.txt
├── protein/      code/ (01–03, structlib.py)   outputs/ (table_p*, fig_p4*, egfr_4hjo_interactive.html)
├── single_cell/  code/ (01–04, common.py)      outputs/ (table_s*, fig_s*, umap_*, dotplot_*)
└── spatial/      code/ (01–04, common.py)      outputs/ (table_t*, fig_t*)
```

The `*.h5ad` files in `outputs/` are intermediate objects (about 150 MB) that the scripts regenerate. Leave them out when you zip the submission.

## Before submitting

Each Agent Audit in `report/report.tex` (sub-questions 1(e), 2(f), 3(f)) ends with a highlighted **"Student verification"** box. The rubric asks what *you* independently verified, so do those checks yourself and write your own result into each box. Then run `bash report/build_report.sh` to rebuild the PDF.
