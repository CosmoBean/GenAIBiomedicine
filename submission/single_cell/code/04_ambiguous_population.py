"""S5: the naive / memory CD4 T boundary (baseline clusters 0 and 1).

Questions:
  1. Are clusters 0 and 1 separated by a gap or by a continuum?
  2. Which markers differ, and how strongly?
  3. Does the 'Naive CD4 T' cluster also contain naive CD8 T cells?
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from common import OUT, MARKERS

sc.settings.figdir = OUT
b = sc.read_h5ad(OUT / "pbmc3k_baseline.h5ad")
raw = b.raw.to_adata()
T = b[b.obs.leiden.isin(["0", "1", "5"])].copy()          # memory CD4, naive CD4, CD8 T
Traw = raw[T.obs_names]


def expr(ad, g):
    x = ad[:, g].X
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x).ravel()


# --- 1. direct DE between the two CD4 clusters
sc.tl.rank_genes_groups(b, "leiden", groups=["1"], reference="0", method="wilcoxon",
                        use_raw=True, key_added="n_vs_m", pts=True)
de = sc.get.rank_genes_groups_df(b, group="1", key="n_vs_m")
de.to_csv(OUT / "table_s5a_DE_cluster1_vs_cluster0.csv", index=False)
show = pd.concat([de.head(12), de.tail(12)])[["names", "logfoldchanges", "pvals_adj", "pct_nz_group"]]
print("cluster 1 (naive) vs 0 (memory): top up / top down\n", show.round(3).to_string())

# --- 2. detection rates of key markers in clusters 0, 1, 5
genes = ["CD3E", "CD4", "CD8A", "CD8B", "CCR7", "SELL", "LEF1", "TCF7", "S100A4", "IL32",
         "AQP3", "CD40LG", "GZMK", "CCL5", "NKG7"]
det = pd.DataFrame({g: expr(Traw, g) > 0 for g in genes}, index=T.obs_names)
det_tbl = det.groupby(T.obs.leiden.astype(str).map(
    {"0": "0 Memory CD4 T", "1": "1 Naive CD4 T", "5": "5 CD8 T"})).mean().T.round(3)
det_tbl.to_csv(OUT / "table_s5b_marker_detection_T_clusters.csv", index_label="gene")
print(det_tbl)

# --- 3. continuum vs gap: naive-minus-memory score along the CD4 clusters
cd4 = T[T.obs.leiden.isin(["0", "1"])].copy()
sc.tl.score_genes(cd4, MARKERS["Naive CD4 T"], score_name="naive", use_raw=True, random_state=0)
sc.tl.score_genes(cd4, MARKERS["Memory CD4 T"], score_name="memory", use_raw=True, random_state=0)
cd4.obs["naive_minus_memory"] = cd4.obs.naive - cd4.obs.memory
d = cd4.obs.naive_minus_memory
# bimodality coefficient (Sarle): > 0.555 suggests bimodality
from scipy.stats import skew, kurtosis
n = len(d); g1 = skew(d); g2 = kurtosis(d)
bc = (g1 ** 2 + 1) / (g2 + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))

# Neighbourhood mixing: fraction of each cell's kNN that is in the other CD4 cluster
conn = b.obsp["connectivities"]
idx = {c: i for i, c in enumerate(b.obs_names)}
lab = b.obs.leiden.astype(str).values
mix = {}
for c in ["0", "1"]:
    other = "1" if c == "0" else "0"
    rows = np.where(lab == c)[0]
    sub = conn[rows]
    tot = np.asarray(sub.sum(1)).ravel()
    oth = np.asarray(sub[:, lab == other].sum(1)).ravel()
    mix[c] = float(np.mean(oth / np.maximum(tot, 1e-9)))

# Naive CD8 contamination: CD8A/CD8B+ cells inside cluster 1 that are CCR7+ and GZMK/CCL5-negative
c1 = Traw[T.obs.leiden == "1"]
cd8p = (expr(c1, "CD8A") > 0) | (expr(c1, "CD8B") > 0)
naive_cd8 = cd8p & (expr(c1, "CCR7") > 0) & (expr(c1, "GZMK") == 0) & (expr(c1, "CCL5") == 0)
cd4p = expr(c1, "CD4") > 0
summary = pd.Series(dict(
    n_cluster0=int((T.obs.leiden == "0").sum()), n_cluster1=int((T.obs.leiden == "1").sum()),
    bimodality_coefficient=bc,
    knn_weight_to_other_cluster_from0=mix["0"], knn_weight_to_other_cluster_from1=mix["1"],
    cluster1_CD8A_or_CD8B_pos=float(cd8p.mean()),
    cluster1_CD8pos_CCR7pos_GZMKneg_CCL5neg=float(naive_cd8.mean()),
    cluster1_CD4_detected=float(cd4p.mean())))
# How well does the score separate the two clusters? (BC alone cannot rule out two close modes)
from sklearn.metrics import roc_auc_score
summary["auc_score_separates_cluster1_from_0"] = roc_auc_score(cd4.obs.leiden == "1", d)
# CD8A|CD8B detection in every baseline cluster, to compare cluster 1 with background (non-T)
# clusters and with true CD8 T (cluster 5). Detection-adjusted CD8 share of cluster 1 =
# (rate in 1 - background) / (rate in 5 - background).
cd8_det = pd.Series({cl: float(((expr(raw[b.obs.leiden == cl], "CD8A") > 0)
                                | (expr(raw[b.obs.leiden == cl], "CD8B") > 0)).mean())
                     for cl in b.obs.leiden.cat.categories}, name="CD8A_or_CD8B_detected")
cd8_det.to_frame().assign(label=b.obs.groupby("leiden", observed=True).cell_type.first().astype(str)).round(3) \
    .to_csv(OUT / "table_s5d_CD8_detection_all_clusters.csv", index_label="cluster")
for bg_name, bg in [("B", cd8_det["2"]), ("NK", cd8_det["4"])]:
    summary[f"cluster1_CD8_share_adjusted_bg_{bg_name}"] = (cd8_det["1"] - bg) / (cd8_det["5"] - bg)
summary.round(3).to_csv(OUT / "table_s5c_ambiguity_summary.csv", header=["value"])
print(summary.round(3))

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for c, col in [("0", "#1f77b4"), ("1", "#ff7f0e")]:
    axes[0].hist(d[cd4.obs.leiden == c], bins=50, alpha=0.6, color=col,
                 label={"0": "cluster 0 (Memory CD4 T)", "1": "cluster 1 (Naive CD4 T)"}[c])
axes[0].set_xlabel("naive score - memory score"); axes[0].set_ylabel("cells")
axes[0].set_title(f"Overlapping distributions\nbimodality coef. = {bc:.2f} (>0.555 = bimodal)")
axes[0].legend(fontsize=8)
order = ["CCR7", "SELL", "LEF1", "TCF7", "S100A4", "IL32", "AQP3", "CD8B"]
dat = [[expr(Traw[T.obs.leiden == c], g) for c in ["1", "0"]] for g in order]
pos = np.arange(len(order))
for j, (c, col) in enumerate([("1", "#ff7f0e"), ("0", "#1f77b4")]):
    vp = axes[1].violinplot([x[j] for x in dat], positions=pos + (j - 0.5) * 0.35, widths=0.33,
                            showmedians=True)
    for body in vp["bodies"]:
        body.set_facecolor(col); body.set_alpha(0.6)
    for part in ("cbars", "cmins", "cmaxes", "cmedians"):
        vp[part].set_color(col)
axes[1].set_xticks(pos, order, rotation=45); axes[1].set_ylabel("log-normalized expression")
axes[1].set_title("Orange = cluster 1 (naive), blue = cluster 0 (memory)")
T.obs["naive_minus_memory"] = np.nan
T.obs.loc[cd4.obs_names, "naive_minus_memory"] = d
um = T.obsm["X_umap"]
m = T.obs.leiden.isin(["0", "1"]).values
lo, hi = np.percentile(d.values, [5, 95])            # centre the colour range on the data
s = axes[2].scatter(um[m, 0], um[m, 1], c=d.values, cmap="viridis", s=6, vmin=lo, vmax=hi)
axes[2].set_title("Clusters 0+1 on UMAP: a gradient (yellow = more naive)\nwith no gap between clusters"); axes[2].set_xticks([]); axes[2].set_yticks([])
plt.colorbar(s, ax=axes[2], label="naive - memory")
fig.tight_layout(); fig.savefig(OUT / "fig_s5a_naive_memory_continuum.png", dpi=130); plt.close(fig)

sc.pl.umap(b, color=["CCR7", "S100A4", "CD8B", "CD4"], use_raw=True, cmap="Reds", ncols=4,
           save="_s5b_T_marker_umaps.png", show=False)
