"""S3: Leiden resolution sweep on the *same* kNN graph (all other settings fixed)."""
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

from common import OUT, RESOLUTIONS, BASELINE_RES, SEED, auto_annotate

sc.settings.figdir = OUT
b = sc.read_h5ad(OUT / "pbmc3k_baseline.h5ad")
keys = [f"leiden_{r}" for r in RESOLUTIONS]

# 1. cluster counts + auto labels at each resolution
rows, labelled = [], {}
for r, k in zip(RESOLUTIONS, keys):
    labels, _ = auto_annotate(b, k)
    b.obs[f"label_{r}"] = b.obs[k].map(labels).astype(str)
    labelled[r] = labels
    for cl, lab in labels.items():
        rows.append(dict(resolution=r, cluster=cl, n_cells=int((b.obs[k] == cl).sum()), label=lab))
lab_tbl = pd.DataFrame(rows)
lab_tbl.to_csv(OUT / "table_s3a_clusters_and_labels_by_resolution.csv", index=False)
print(lab_tbl.groupby("resolution").size().rename("n_clusters"))

# 2. pairwise ARI / NMI between partitions
M = pd.DataFrame(index=RESOLUTIONS, columns=RESOLUTIONS, dtype=float)
N = M.copy()
for r1, r2 in itertools.product(RESOLUTIONS, RESOLUTIONS):
    M.loc[r1, r2] = ARI(b.obs[f"leiden_{r1}"], b.obs[f"leiden_{r2}"])
    N.loc[r1, r2] = NMI(b.obs[f"leiden_{r1}"], b.obs[f"leiden_{r2}"])
M.round(3).to_csv(OUT / "table_s3b_pairwise_ARI.csv"); N.round(3).to_csv(OUT / "table_s3c_pairwise_NMI.csv")
print("ARI\n", M.round(3), "\nNMI\n", N.round(3))

# 3. splits / merges between consecutive resolutions (>=10% of a parent cluster)
events = []
for lo, hi in zip(RESOLUTIONS[:-1], RESOLUTIONS[1:]):
    ct = pd.crosstab(b.obs[f"leiden_{lo}"], b.obs[f"leiden_{hi}"])
    frac = ct.div(ct.sum(1), axis=0)
    for parent, row in frac.iterrows():
        kids = row[row >= 0.10].index.tolist()
        if len(kids) > 1:
            events.append(dict(
                from_res=lo, to_res=hi, parent=f"{parent} ({labelled[lo][parent]}, n={ct.loc[parent].sum()})",
                children="; ".join(f"{c} ({labelled[hi][c]}, {row[c]:.0%})" for c in kids)))
ev = pd.DataFrame(events)
ev.to_csv(OUT / "table_s3d_splits_between_resolutions.csv", index=False)
print(ev.to_string())

# 4. contingency heatmaps vs baseline
fig, axes = plt.subplots(1, len(RESOLUTIONS) - 1, figsize=(5.2 * (len(RESOLUTIONS) - 1), 4.6))
others = [r for r in RESOLUTIONS if r != BASELINE_RES]
for ax, r in zip(axes, others):
    ct = pd.crosstab(b.obs[f"label_{BASELINE_RES}"], b.obs[f"leiden_{r}"])
    sns.heatmap(ct.div(ct.sum(1), axis=0), cmap="Blues", vmin=0, vmax=1, annot=ct, fmt="d",
                annot_kws={"size": 7}, cbar=False, ax=ax)
    ax.set_title(f"baseline label vs leiden res={r}\nARI={M.loc[BASELINE_RES, r]:.2f}  NMI={N.loc[BASELINE_RES, r]:.2f}")
    ax.set_xlabel(f"cluster @ res {r}"); ax.set_ylabel("baseline (res 1.0) label")
fig.tight_layout(); fig.savefig(OUT / "fig_s3a_contingency_vs_baseline.png", dpi=130); plt.close(fig)

sc.pl.umap(b, color=keys, legend_loc="on data", ncols=4, save="_s3b_resolutions.png", show=False)

# 5. per-population label stability: fraction of baseline-labelled cells that
#    receive the *same* label at each other resolution
stab = {}
for r in RESOLUTIONS:
    same = pd.Series(b.obs[f"label_{r}"].astype(str).to_numpy() == b.obs[f"label_{BASELINE_RES}"].astype(str).to_numpy(), index=b.obs_names)
    stab[r] = same.groupby(b.obs[f"label_{BASELINE_RES}"].astype(str)).mean()
stab = pd.DataFrame(stab).round(3)
stab.to_csv(OUT / "table_s3e_label_stability_by_resolution.csv", index_label="baseline_label")
print(stab)

# 6. extra check: Leiden seed sensitivity at the baseline resolution
# Also: per-label fraction of baseline cells keeping their label, and one run with seed 0
# iterated to convergence (n_iterations=-1) to separate seed noise from early stopping.
seed_rows = []
base_lab = b.obs.cell_type.astype(str)
for s, n_it in [(1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (SEED, -1)]:
    sc.tl.leiden(b, resolution=BASELINE_RES, key_added="tmp", random_state=s,
                 flavor="igraph", n_iterations=n_it, directed=False)
    lab_s, _ = auto_annotate(b, "tmp")
    kept = (b.obs.tmp.map(lab_s).astype(str) == base_lab).groupby(base_lab).mean()
    seed_rows.append(dict(seed=s, n_iterations=n_it, n_clusters=b.obs.tmp.nunique(),
                          ARI_vs_seed0=ARI(b.obs.tmp, b.obs[f"leiden_{BASELINE_RES}"]),
                          **{f"kept_{k}": v for k, v in kept.round(3).items()}))
seeds = pd.DataFrame(seed_rows)
seeds.to_csv(OUT / "table_s3f_leiden_seed_sensitivity.csv", index=False)
print(seeds)
