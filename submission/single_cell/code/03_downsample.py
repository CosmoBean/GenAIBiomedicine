"""S4: fixed-seed 70% cell subsample, full pipeline rerun, compare on shared cells."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

from common import (OUT, BASELINE_RES, DOWNSAMPLE_SEED, DOWNSAMPLE_FRAC,
                    load_raw, apply_qc, preprocess_and_cluster, auto_annotate)

sc.settings.figdir = OUT
base = sc.read_h5ad(OUT / "pbmc3k_baseline.h5ad")
base_lab = base.obs["cell_type"].astype(str)

# subsample QC-passed cells (same QC as baseline), then rerun normalization,
# HVG selection, PCA, kNN, UMAP and Leiden from raw counts
a = apply_qc(load_raw())
rng = np.random.default_rng(DOWNSAMPLE_SEED)
keep = rng.choice(a.obs_names, size=int(round(DOWNSAMPLE_FRAC * a.n_obs)), replace=False)
sub = preprocess_and_cluster(a[keep].copy(), resolutions=(BASELINE_RES,))
sub.obs["leiden"] = sub.obs[f"leiden_{BASELINE_RES}"]
labels, scores = auto_annotate(sub, "leiden")
sub.obs["cell_type"] = sub.obs["leiden"].map(labels).astype("category")
print(f"kept {sub.n_obs}/{a.n_obs} cells (seed={DOWNSAMPLE_SEED}); "
      f"{sub.obs.leiden.nunique()} clusters; HVG overlap with baseline: "
      f"{len(set(sub.var_names) & set(base.var_names))}/{sub.n_vars}")

shared = sub.obs_names.intersection(base.obs_names)
bc, sc_ = base.obs.loc[shared, "leiden"], sub.obs.loc[shared, "leiden"]
bl, sl = base_lab.loc[shared], sub.obs.loc[shared, "cell_type"].astype(str)
summary = pd.Series(dict(
    n_baseline=base.n_obs, n_subsample=sub.n_obs, n_shared=len(shared),
    n_clusters_baseline=base.obs.leiden.nunique(), n_clusters_subsample=sub.obs.leiden.nunique(),
    ARI_clusters=ARI(bc, sc_), NMI_clusters=NMI(bc, sc_),
    ARI_labels=ARI(bl, sl), label_agreement=(bl == sl).mean()))
summary.to_csv(OUT / "table_s4a_downsample_summary.csv", header=["value"])
print(summary)

per = (bl == sl).groupby(bl).agg(["mean", "size"]).rename(columns={"mean": "same_label_frac", "size": "n_shared"})
per.round(3).to_csv(OUT / "table_s4b_downsample_label_stability.csv", index_label="baseline_label")
print(per.round(3))

ct = pd.crosstab(bl, sl)
fig, ax = plt.subplots(figsize=(7.5, 5.5))
sns.heatmap(ct.div(ct.sum(1), axis=0), annot=ct, fmt="d", cmap="Blues", vmin=0, vmax=1, cbar=False, ax=ax)
ax.set_xlabel(f"subsample label ({DOWNSAMPLE_FRAC:.0%}, seed {DOWNSAMPLE_SEED})"); ax.set_ylabel("baseline label")
ax.set_title(f"Shared cells (n={len(shared)}): ARI={summary.ARI_labels:.2f}, agreement={summary.label_agreement:.1%}")
fig.tight_layout(); fig.savefig(OUT / "fig_s4a_downsample_confusion.png", dpi=130); plt.close(fig)

# marker reproducibility: top-20 Wilcoxon markers per label, Jaccard between runs
def top_markers(ad, n=20):
    ad.obs["_lab"] = ad.obs["cell_type"].astype(str)
    sc.tl.rank_genes_groups(ad, "_lab", method="wilcoxon", use_raw=True, key_added="_rg")
    df = sc.get.rank_genes_groups_df(ad, group=None, key="_rg")
    return {g: set(d.sort_values("scores", ascending=False).names[:n]) for g, d in df.groupby("group")}

mb, ms = top_markers(base), top_markers(sub)
jac = pd.Series({k: len(mb[k] & ms[k]) / len(mb[k] | ms[k]) for k in mb if k in ms}, name="jaccard_top20")
jac.round(3).to_csv(OUT / "table_s4c_marker_jaccard.csv", index_label="label")
print(jac.round(2))

sub.obs["baseline_label"] = base_lab.reindex(sub.obs_names).values
sc.pl.umap(sub, color=["leiden", "cell_type", "baseline_label"], ncols=3, wspace=0.5,
           save="_s4b_downsample.png", show=False)
sub.write_h5ad(OUT / "pbmc3k_downsample70.h5ad", compression="gzip")

# --- additional checks (not part of the primary seed-42 comparison) ---
# (a) Is the lost CD8 T / NK split just a resolution-scaling effect? Sweep
#     resolution on the same seed-42 subsample.
from common import RESOLUTIONS
rows = []
for r in RESOLUTIONS + [1.2, 1.6, 2.0]:
    sc.tl.leiden(sub, resolution=r, key_added="tmp", random_state=0, flavor="igraph",
                 n_iterations=2, directed=False)
    lab, _ = auto_annotate(sub, "tmp")
    l = sub.obs.tmp.map(lab).astype(str)
    rows.append(dict(resolution=r, n_clusters=sub.obs.tmp.nunique(),
                     ARI_vs_baseline=ARI(base.obs.loc[shared, "leiden"], sub.obs.loc[shared, "tmp"]),
                     label_agreement=(l.loc[shared] == bl).mean(),
                     labels="; ".join(sorted(set(lab.values())))))
rs = pd.DataFrame(rows).sort_values("resolution")
rs.to_csv(OUT / "table_s4d_subsample_resolution_sweep.csv", index=False)
print(rs.to_string())

# (b) Other subsample seeds at the baseline resolution
rows = []
for s in [1, 2, 3, 4, 5]:
    k = np.random.default_rng(s).choice(a.obs_names, size=len(keep), replace=False)
    t = preprocess_and_cluster(a[k].copy(), resolutions=(BASELINE_RES,))
    t.obs["leiden"] = t.obs[f"leiden_{BASELINE_RES}"]
    lab, _ = auto_annotate(t, "leiden")
    l = t.obs.leiden.map(lab).astype(str)
    sh = t.obs_names
    rows.append(dict(seed=s, n_clusters=t.obs.leiden.nunique(),
                     ARI_vs_baseline=ARI(base.obs.loc[sh, "leiden"], t.obs.leiden),
                     label_agreement=(l == base_lab.loc[sh]).mean(),
                     CD8T_cluster_found="CD8 T" in set(lab.values()),
                     platelet_cluster_found="Platelet" in set(lab.values())))
ms_ = pd.DataFrame(rows)
ms_.to_csv(OUT / "table_s4e_other_subsample_seeds.csv", index=False)
print(ms_.to_string())
