"""S1 + S2: QC, baseline clustering, data-derived markers, annotation."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from common import (OUT, QC, SEED, N_HVG, N_PCS, N_NEIGHBORS, BASELINE_RES,
                    RESOLUTIONS, load_raw, apply_qc, preprocess_and_cluster,
                    auto_annotate)

sc.settings.figdir = OUT
sc.settings.set_figure_params(dpi=120, frameon=False)

# ---------------- QC ----------------
raw = load_raw()
before = dict(cells=raw.n_obs, genes=raw.n_vars)

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
for ax, col, thr in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                        [(QC["min_genes"], QC["max_genes"]), None, (QC["max_pct_mt"],)]):
    ax.hist(raw.obs[col], bins=80, color="#4c72b0")
    for t in (thr or []):
        ax.axvline(t, color="crimson", ls="--")
    med = raw.obs[col].median()
    mad = 1.4826 * np.median(np.abs(raw.obs[col] - med))
    for k in (3, 5):
        ax.axvline(med + k * mad, color="grey", ls=":", lw=1)
    ax.set_title(col)
    ax.set_ylabel("cells")
fig.suptitle("PBMC3k QC (red dashed = chosen thresholds; grey dotted = median + 3/5 MAD)")
fig.tight_layout()
fig.savefig(OUT / "fig_s1a_qc_histograms.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
sc_ = ax.scatter(raw.obs.total_counts, raw.obs.n_genes_by_counts, c=raw.obs.pct_counts_mt,
                 s=4, cmap="viridis", vmax=10)
ax.axhline(QC["max_genes"], color="crimson", ls="--")
ax.set_xlabel("total UMI counts"); ax.set_ylabel("genes detected")
plt.colorbar(sc_, label="% mitochondrial")
fig.tight_layout(); fig.savefig(OUT / "fig_s1b_qc_scatter.png"); plt.close(fig)

a = apply_qc(raw)
after = dict(cells=a.n_obs, genes=a.n_vars)
removed = dict(
    low_genes=int((raw.obs.n_genes_by_counts < QC["min_genes"]).sum()),
    high_genes=int((raw.obs.n_genes_by_counts >= QC["max_genes"]).sum()),
    high_mt=int((raw.obs.pct_counts_mt >= QC["max_pct_mt"]).sum()),
)
print("QC before", before, "after", after, "removed (overlapping)", removed)

# ---------------- baseline ----------------
b = preprocess_and_cluster(a, resolutions=RESOLUTIONS)
b.obs["leiden"] = b.obs[f"leiden_{BASELINE_RES}"]

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.plot(np.arange(1, 51), b.uns["pca"]["variance_ratio"], "o-", ms=3)
ax.axvline(N_PCS, color="crimson", ls="--", label=f"n_pcs = {N_PCS}")
ax.set_xlabel("PC"); ax.set_ylabel("variance ratio"); ax.set_yscale("log"); ax.legend()
fig.tight_layout(); fig.savefig(OUT / "fig_s1c_pca_elbow.png"); plt.close(fig)

sizes = b.obs["leiden"].value_counts().sort_index(key=lambda s: s.astype(int))
sizes.rename("n_cells").to_csv(OUT / "table_s1_baseline_cluster_sizes.csv", index_label="cluster")
print("baseline clusters:", len(sizes)); print(sizes.to_string())

sc.pl.umap(b, color="leiden", legend_loc="on data", title=f"Baseline Leiden (res={BASELINE_RES})",
           save="_s1d_baseline_clusters.png", show=False)

# ---------------- data-derived markers ----------------
sc.tl.rank_genes_groups(b, "leiden", method="wilcoxon", use_raw=True, pts=True)
rg = sc.get.rank_genes_groups_df(b, group=None)
rg.to_csv(OUT / "table_s2a_wilcoxon_all_markers.csv", index=False)
top = (rg[(rg.pvals_adj < 0.05) & (rg.logfoldchanges > 1)]
       .sort_values(["group", "scores"], ascending=[True, False])
       .groupby("group").head(10))
top_tbl = top.groupby("group").apply(
    lambda d: ", ".join(f"{g} ({lfc:.1f})" for g, lfc in zip(d.names, d.logfoldchanges)),
    include_groups=False)
top_tbl = top_tbl.reindex(sorted(top_tbl.index, key=int))
top_tbl.rename("top10 data-derived markers (log2FC)").to_csv(OUT / "table_s2b_top10_markers.csv",
                                                              index_label="cluster")
print(top_tbl.to_string())

sc.pl.rank_genes_groups_dotplot(b, n_genes=4, groupby="leiden", standard_scale="var",
                                save="_s2a_top_markers_by_cluster.png", show=False)

# ---------------- annotation ----------------
labels, scores = auto_annotate(b, "leiden")
scores.round(2).to_csv(OUT / "table_s2c_marker_set_scores_baseline.csv", index_label="cluster")
b.obs["cell_type"] = b.obs["leiden"].map(labels).astype("category")
ann = pd.DataFrame({"cluster": list(labels), "label": list(labels.values()),
                    "n_cells": [int(sizes[c]) for c in labels]})
ann = ann.merge(top_tbl.rename("top_markers").reset_index().rename(columns={"group": "cluster"}),
                on="cluster", how="left")
ann.to_csv(OUT / "table_s2d_baseline_annotation.csv", index=False)
print(ann.to_string())

sc.pl.umap(b, color="cell_type", title="Baseline annotation (marker-score rule)",
           save="_s2b_baseline_annotated.png", show=False)

from common import MARKERS
genes = [g for v in MARKERS.values() for g in v[:3] if g in b.raw.var_names]
sc.pl.dotplot(b, list(dict.fromkeys(genes)), groupby="leiden", standard_scale="var",
              save="_s2c_canonical_markers.png", show=False)
sc.pl.umap(b, color=["IL7R", "CCR7", "S100A4", "CD8A", "NKG7", "MS4A1", "CD14",
                     "FCGR3A", "FCER1A", "PPBP"], ncols=5, cmap="Reds", use_raw=True,
           save="_s2d_marker_feature_plots.png", show=False)

b.write_h5ad(OUT / "pbmc3k_baseline.h5ad", compression="gzip")
json.dump(dict(seed=SEED, qc=QC, before=before, after=after, removed=removed,
               normalize="normalize_total(1e4)+log1p", hvg=f"seurat top {N_HVG}",
               regress_out=["total_counts", "pct_counts_mt"], scale_max=10,
               n_pcs=N_PCS, n_neighbors=N_NEIGHBORS, leiden="igraph, n_iterations=2",
               baseline_resolution=BASELINE_RES, n_clusters=int(len(sizes)),
               scanpy=sc.__version__),
          open(OUT / "baseline_settings.json", "w"), indent=2)
