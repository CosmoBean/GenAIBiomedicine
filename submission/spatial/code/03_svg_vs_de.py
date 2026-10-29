"""T3 + T4: spatially variable genes (Moran's I) vs differential expression.

Gene set (same for both rankings): the top N_SVG_GENES genes by seurat_v3
HVG rank computed in 01 from expression only (mitochondrial genes excluded).

Spatial ranking : Moran's I on the 6-NN hex-grid graph (row-normalised
                  weights), analytic p-value, Benjamini-Hochberg FDR.
DE ranking      : Wilcoxon one-vs-rest across expression domains; each
                  gene's DE strength = its best (max) Wilcoxon z-score over
                  domains. BH FDR from scanpy. Because the domains were defined from
                  these same HVGs, DE p-values are circular and are not used as evidence.
Sensitivity     : the DE ranking is repeated with effect size (max log2FC over
                  domains) instead of the Wilcoxon z, because z is rank-based,
                  grows with group size and is penalised by zeros inside a domain.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import spearmanr

from common import OUT, SEED, N_SVG_GENES

a = sc.read_h5ad(OUT / "visium_domains.h5ad")   # log-normalised X, grid graph from 02
hv = a.var[a.var.highly_variable & ~a.var.mt].sort_values("hvg_rank")
genes = hv.index[:N_SVG_GENES].tolist()
pd.Series(genes, name="gene").to_csv(OUT / "table_t3a_tested_gene_set.csv", index=False)

# ---------- spatial ranking ----------
sq.gr.spatial_autocorr(a, mode="moran", genes=genes, connectivity_key="spatial_connectivities",
                       transformation=True, n_perms=None, seed=SEED, show_progress_bar=False)
mor = a.uns["moranI"].loc[genes, ["I", "pval_norm", "pval_norm_fdr_bh"]].copy()
mor["spatial_rank"] = mor.I.rank(ascending=False).astype(int)

# ---------- DE ranking (same genes) ----------
sc.tl.rank_genes_groups(a, "domain", method="wilcoxon", use_raw=False, key_added="de")
de = sc.get.rank_genes_groups_df(a, group=None, key="de")
de = de[de.names.isin(genes)]
best = de.loc[de.groupby("names").scores.idxmax()].set_index("names")
mor["DE_best_z"] = best.scores
mor["DE_best_domain"] = best.group
mor["DE_best_log2FC"] = best.logfoldchanges
mor["DE_padj"] = best.pvals_adj
mor["DE_rank"] = mor.DE_best_z.rank(ascending=False).astype(int)
mor["DE_max_log2FC"] = de.groupby("names").logfoldchanges.max()
mor["DE_rank_log2FC"] = mor.DE_max_log2FC.rank(ascending=False, method="min").astype(int)
mor["mean_log_expr"] = np.asarray(a[:, genes].X.mean(0)).ravel()
mor["frac_spots_detected"] = np.asarray((a[:, genes].X > 0).mean(0)).ravel()

rho, p = spearmanr(mor.I, mor.DE_best_z)
rho_fc, p_fc = spearmanr(mor.I, mor.DE_max_log2FC)
pd.DataFrame([dict(DE_metric="max Wilcoxon z over domains (main)", spearman_rho=rho, p=p),
              dict(DE_metric="max log2FC over domains (sensitivity)", spearman_rho=rho_fc, p=p_fc)]
             ).round(4).to_csv(OUT / "table_t4b_rank_agreement_by_DE_metric.csv", index=False)
top = 30
both = mor[(mor.spatial_rank <= top) & (mor.DE_rank <= top)].sort_values("spatial_rank")
de_only = mor[(mor.DE_rank <= top) & (mor.spatial_rank > 2 * top)].sort_values("DE_rank")
sp_only = mor[(mor.spatial_rank <= top) & (mor.DE_rank > 2 * top)].sort_values("spatial_rank")
mor["category"] = "other"
mor.loc[both.index, "category"] = f"top{top} in both"
mor.loc[de_only.index, "category"] = f"DE top{top}, spatial rank >{2*top}"
mor.loc[sp_only.index, "category"] = f"spatial top{top}, DE rank >{2*top}"
mor.sort_values("spatial_rank").round(4).to_csv(OUT / "table_t3b_t4_moran_and_DE_rankings.csv", index_label="gene")

print(f"{len(genes)} genes; Moran's I FDR<0.05: {(mor.pval_norm_fdr_bh < 0.05).sum()}; "
      f"Spearman(I, DE z) = {rho:.2f} (p={p:.1e}); Spearman(I, max log2FC) = {rho_fc:.2f}")
print("top 15 spatial:\n", mor.sort_values("spatial_rank").head(15)[["I", "pval_norm_fdr_bh", "DE_rank", "DE_best_domain"]].round(3))
for name, d in [("BOTH", both), ("DE-high / spatial-low", de_only), ("spatial-high / DE-low", sp_only)]:
    print(name, "\n", d[["I", "spatial_rank", "DE_best_z", "DE_rank", "DE_best_domain", "DE_best_log2FC",
                         "frac_spots_detected"]].round(3).head(8))

# ---------- figures ----------
fig, ax = plt.subplots(figsize=(6.5, 5.5))
colors = {f"top{top} in both": "crimson", f"DE top{top}, spatial rank >{2*top}": "royalblue",
          f"spatial top{top}, DE rank >{2*top}": "darkorange", "other": "#bbbbbb"}
for cat, col in colors.items():
    d = mor[mor.category == cat]
    ax.scatter(d.DE_best_z, d.I, s=14, c=col, label=f"{cat} (n={len(d)})")
lab = pd.concat([both.head(4), de_only.head(4), sp_only.head(4)])
for g, row in lab.iterrows():
    ax.annotate(g, (row.DE_best_z, row.I), fontsize=7, xytext=(3, 2), textcoords="offset points")
ax.set_xlabel("DE strength: max one-vs-rest Wilcoxon z over domains")
ax.set_ylabel("Moran's I (6-NN hex grid)")
ax.set_title(f"Same {len(genes)} HVGs ranked two ways (Spearman ρ = {rho:.2f})")
ax.legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(OUT / "fig_t4_DE_vs_moran.png", dpi=140); plt.close(fig)

show = (mor.sort_values("spatial_rank").index[:3].tolist()
        + de_only.index[:2].tolist() + sp_only.index[:1].tolist())
show = list(dict.fromkeys(show))
sq.pl.spatial_scatter(a, color=show, ncols=3, size=1.3, cmap="viridis", img_alpha=0.4,
                      title=[f"{g}  I={mor.loc[g,'I']:.2f}, DE rank {mor.loc[g,'DE_rank']}" for g in show])
plt.savefig(OUT / "fig_t3_selected_genes_on_tissue.png", dpi=120, bbox_inches="tight"); plt.close("all")

# ---------- technical-confound check ----------
# Is sequencing depth itself spatially structured, and does it track domains?
for v in ["total_counts", "pct_counts_mt"]:
    a.obs[v + "_log"] = np.log1p(a.obs[v])
tech = sc.AnnData(a.obs[["total_counts_log", "pct_counts_mt"]].to_numpy(), obs=a.obs[[]].copy())
tech.var_names = ["log_total_counts", "pct_counts_mt"]
tech.obsp["spatial_connectivities"] = a.obsp["spatial_connectivities"]
sq.gr.spatial_autocorr(tech, mode="moran", connectivity_key="spatial_connectivities",
                       n_perms=None, show_progress_bar=False)
eta2 = {}
for v in ["total_counts_log", "pct_counts_mt"]:
    x = a.obs[v]; g = a.obs.domain
    eta2[v] = float(((x.groupby(g).transform("mean") - x.mean()) ** 2).sum() / ((x - x.mean()) ** 2).sum())
tc = tech.uns["moranI"][["I", "pval_norm"]].copy()
tc["eta2_explained_by_domain"] = [eta2["total_counts_log"], eta2["pct_counts_mt"]]
tc.round(4).to_csv(OUT / "table_t5b_technical_covariate_check.csv", index_label="covariate")
print(tc.round(3))
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, v in zip(axes, ["total_counts_log", "pct_counts_mt"]):
    xy = a.obsm["spatial"]
    s = ax.scatter(xy[:, 0], -xy[:, 1], c=a.obs[v], s=5, cmap="magma"); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{v}: Moran's I={tc.iloc[0 if v.startswith('total') else 1].I:.2f}, "
                 f"η² by domain={eta2[v]:.2f}", fontsize=9)
    plt.colorbar(s, ax=ax, shrink=0.7)
fig.tight_layout(); fig.savefig(OUT / "fig_t5b_technical_covariates_on_tissue.png", dpi=120); plt.close(fig)
