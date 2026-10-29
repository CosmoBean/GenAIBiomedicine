"""Shared settings and helpers for the PBMC3k analysis (Task 2).

Every script imports its parameters from here so that the baseline,
resolution sweep and downsampling run use identical settings.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "filtered_gene_bc_matrices" / "hg19"
OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(exist_ok=True)

SEED = 0                 # PCA / neighbors / Leiden / UMAP
DOWNSAMPLE_SEED = 42     # cell subsample in 03_downsample.py
DOWNSAMPLE_FRAC = 0.70

QC = dict(min_genes=200, max_genes=2500, max_pct_mt=5.0, min_cells=3)
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
BASELINE_RES = 1.0
RESOLUTIONS = [0.3, 0.6, 1.0, 1.4]

# Canonical PBMC marker genes, taken from literature (not from any
# reference label set). Used only *after* clustering, to score clusters.
# Sources are listed in the report (Task 2, Table S2).
MARKERS = {
    "Naive CD4 T": ["IL7R", "CCR7", "LEF1", "TCF7", "SELL"],
    "Memory CD4 T": ["IL7R", "S100A4", "AQP3", "CD40LG", "IL32"],
    "CD8 T": ["CD8A", "CD8B", "CCL5", "GZMK"],
    "NK": ["GNLY", "NKG7", "KLRD1", "GZMB", "PRF1", "FCGR3A"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7", "LST1", "IFITM3", "RHOC"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "CD1C"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C", "TCF4"],
    "Platelet": ["PPBP", "PF4", "GNG11", "SPARC"],
}


def load_raw():
    a = sc.read_10x_mtx(DATA, var_names="gene_symbols", cache=False)
    a.var_names_make_unique()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True,
                               percent_top=None, log1p=False)
    return a


def apply_qc(a):
    a = a[(a.obs.n_genes_by_counts >= QC["min_genes"])
          & (a.obs.n_genes_by_counts < QC["max_genes"])
          & (a.obs.pct_counts_mt < QC["max_pct_mt"])].copy()
    sc.pp.filter_genes(a, min_cells=QC["min_cells"])
    return a


def preprocess_and_cluster(a, resolutions=(BASELINE_RES,), seed=SEED):
    """Normalize -> log1p -> HVG -> scale -> PCA -> kNN -> UMAP -> Leiden.

    `a` must contain raw counts in .X. Full log-normalized data is kept
    in .raw for marker testing/plotting. Returns a new AnnData.
    """
    a = a.copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat")
    a = a[:, a.var.highly_variable].copy()
    sc.pp.regress_out(a, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, random_state=seed)
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    for r in resolutions:
        sc.tl.leiden(a, resolution=r, key_added=f"leiden_{r}", random_state=seed,
                     flavor="igraph", n_iterations=2, directed=False)
    return a


def marker_scores(a, groupby):
    """Mean z-scored (across cells) log-expression of each marker set per cluster."""
    raw = a.raw.to_adata()
    rows = {}
    for label, genes in MARKERS.items():
        g = [x for x in genes if x in raw.var_names]
        X = raw[:, g].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        z = (X - X.mean(0)) / (X.std(0) + 1e-9)
        rows[label] = pd.Series(z.mean(1), index=a.obs_names)
    S = pd.DataFrame(rows)
    return S.groupby(a.obs[groupby].values).mean()


def auto_annotate(a, groupby, margin=0.15):
    """Label each cluster by its top marker-set score.

    If the best and second-best scores differ by < `margin`, the cluster
    is reported as 'A/B (ambiguous)' instead of forcing a label.
    """
    S = marker_scores(a, groupby)
    labels = {}
    for cl, row in S.iterrows():
        top = row.sort_values(ascending=False)
        if top.iloc[0] - top.iloc[1] < margin:
            labels[cl] = f"{top.index[0]}/{top.index[1]} (ambiguous)"
        else:
            labels[cl] = top.index[0]
    return labels, S
