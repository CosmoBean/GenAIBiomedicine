"""T1: expression-only clustering of Visium spots, then projection onto tissue.

Spatial coordinates and the image are *removed* from the object used for
feature selection, PCA, the kNN graph and Leiden. They are re-attached only
after the labels exist. An assertion enforces this.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from common import DATA, COUNTS, OUT, SEED, QC, N_HVG, N_PCS, N_NEIGHBORS, RES

sc.settings.figdir = OUT
sc.settings.set_figure_params(dpi=120, frameon=False)

full = sq.read.visium(DATA, counts_file=COUNTS)
full.var_names_make_unique()
full.var["mt"] = full.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(full, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
before = dict(spots=full.n_obs, genes=full.n_vars, in_tissue=int(full.obs.in_tissue.sum()))

fig, axes = plt.subplots(1, 3, figsize=(13, 3.4))
for ax, c, t in zip(axes, ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
                    [QC["min_counts"], QC["min_genes"], QC["max_pct_mt"]]):
    ax.hist(full.obs[c], bins=80); ax.axvline(t, color="crimson", ls="--"); ax.set_title(c)
fig.tight_layout(); fig.savefig(OUT / "fig_t0_qc_histograms.png"); plt.close(fig)

keep = ((full.obs.total_counts >= QC["min_counts"]) & (full.obs.n_genes_by_counts >= QC["min_genes"])
        & (full.obs.pct_counts_mt < QC["max_pct_mt"]) & (full.obs.in_tissue == 1))
full = full[keep].copy()
sc.pp.filter_genes(full, min_cells=QC["min_cells"])
after = dict(spots=full.n_obs, genes=full.n_vars)
print("QC", before, "->", after)

# -------- expression-only object: strip every spatial field --------
expr = sc.AnnData(X=full.X.copy(), obs=full.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy(),
                  var=full.var[["gene_ids", "mt"]].copy())
assert "spatial" not in expr.obsm and "spatial" not in expr.uns
assert not {"array_row", "array_col", "in_tissue"} & set(expr.obs.columns)

expr.layers["counts"] = expr.X.copy()
sc.pp.highly_variable_genes(expr, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts")
sc.pp.normalize_total(expr, target_sum=1e4)
sc.pp.log1p(expr)
expr.raw = expr
hv = expr[:, expr.var.highly_variable].copy()
sc.pp.scale(hv, max_value=10)
sc.tl.pca(hv, n_comps=50, random_state=SEED)
sc.pp.neighbors(hv, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)
sc.tl.umap(hv, random_state=SEED)
sc.tl.leiden(hv, resolution=RES, key_added="domain", random_state=SEED,
             flavor="igraph", n_iterations=2, directed=False)

# -------- re-attach coordinates only now --------
full.obs["domain"] = hv.obs["domain"].reindex(full.obs_names)
full.obsm["X_umap"] = hv.obsm["X_umap"]
full.obsm["X_pca"] = hv.obsm["X_pca"]
full.layers["counts"] = full.X.copy()
sc.pp.normalize_total(full, target_sum=1e4); sc.pp.log1p(full)
full.var["highly_variable"] = expr.var["highly_variable"]
full.var["hvg_rank"] = expr.var["highly_variable_rank"]
full.uns["domain_colors"] = None
sizes = full.obs.domain.value_counts().sort_index(key=lambda s: s.astype(int))
sizes.rename("n_spots").to_csv(OUT / "table_t1_domain_sizes.csv", index_label="domain")
print(len(sizes), "domains"); print(sizes.to_string())

# same palette in both spaces
pal = dict(zip(sizes.index, sc.pl.palettes.default_20[:len(sizes)]))
full.uns["domain_colors"] = [pal[d] for d in full.obs.domain.cat.categories]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
sc.pl.umap(full, color="domain", ax=axes[0], show=False, legend_loc="on data",
           title=f"Expression space (UMAP), Leiden res={RES}")
sq.pl.spatial_scatter(full, color="domain", ax=axes[1], title="Same labels at tissue coordinates",
                      size=1.3, legend_loc="right margin")
fig.tight_layout(); fig.savefig(OUT / "fig_t1_domains_expression_vs_tissue.png", dpi=130); plt.close(fig)

# per-domain small multiples on tissue (helps judge contiguity)
cats = list(full.obs.domain.cat.categories)
nc = 5; nr = int(np.ceil(len(cats) / nc))
fig, axes = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.0 * nr))
xy = full.obsm["spatial"]
for ax, d in zip(axes.ravel(), cats):
    m = (full.obs.domain == d).values
    ax.scatter(xy[~m, 0], -xy[~m, 1], s=1, c="#dddddd"); ax.scatter(xy[m, 0], -xy[m, 1], s=1.5, c=pal[d])
    ax.set_title(f"domain {d} (n={m.sum()})", fontsize=9); ax.set_aspect("equal"); ax.axis("off")
for ax in axes.ravel()[len(cats):]:
    ax.axis("off")
fig.tight_layout(); fig.savefig(OUT / "fig_t1b_domains_small_multiples.png", dpi=120); plt.close(fig)

# markers per domain (used for T4 and for domain descriptions)
sc.tl.rank_genes_groups(full, "domain", method="wilcoxon", use_raw=False)
rg = sc.get.rank_genes_groups_df(full, group=None)
rg.to_csv(OUT / "table_t1b_domain_markers_wilcoxon_allgenes.csv", index=False)
top = rg[rg.logfoldchanges > 1].sort_values(["group", "scores"], ascending=[True, False]).groupby("group").head(6)
top.groupby("group").names.apply(", ".join).to_csv(OUT / "table_t1c_domain_top_markers.csv")
print(top.groupby("group").names.apply(", ".join).to_string())

del full.uns["rank_genes_groups"]
full.write_h5ad(OUT / "visium_domains.h5ad", compression="gzip")
json.dump(dict(seed=SEED, qc=QC, before=before, after=after, hvg=f"seurat_v3 top {N_HVG} on raw counts",
               normalize="normalize_total 1e4 + log1p", scale_max=10, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS,
               leiden=f"igraph, res={RES}, n_iterations=2", n_domains=int(len(sizes)),
               spatial_fields_used_for_clustering="none (asserted)"),
          open(OUT / "domain_settings.json", "w"), indent=2)
