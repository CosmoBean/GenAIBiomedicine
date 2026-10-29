"""T2 + T5: spatial organisation of expression-defined domains and a randomised null.

Spatial graph: Visium hexagonal grid on obsm['spatial'] (spot centres; same edges as the array-coordinate hex lattice), first ring
only (each spot -> up to 6 physically adjacent spots, ~100 um centre-to-centre).
Only in-tissue spots that passed QC are nodes, so off-tissue spots and
QC-removed spots can never be neighbours; edge spots simply have < 6 neighbours.

Statistics on the observed labels:
  S1  global same-domain edge fraction  = (# edges joining two spots of the same
      domain) / (# edges). Higher = more spatially coherent domains.
  S2  per-domain self-neighbour fraction (same quantity restricted to one domain)
  S3  squidpy neighbourhood-enrichment z-scores (domain x domain)
  S4  Moran's I of each domain's 0/1 membership indicator

Null A (main, T5): permute domain labels across spots on the FIXED graph
(preserves domain sizes and the graph exactly; destroys spatial placement).
Null B (check): permute spot coordinates, rebuild the grid graph from the
permuted coordinates every repetition (the "shuffle coordinates" control).
Both use N_PERM repetitions with a recorded seed.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import squidpy as sq

from common import OUT, SEED, N_PERM

a = sc.read_h5ad(OUT / "visium_domains.h5ad")
sq.gr.spatial_neighbors(a, coord_type="grid", n_neighs=6, n_rings=1, spatial_key="spatial")
A = a.obsp["spatial_connectivities"].tocsr()
A = ((A + A.T) > 0).astype(np.int8)
deg = np.asarray(A.sum(1)).ravel()
graph_info = dict(nodes=int(a.n_obs), undirected_edges=int(A.sum() // 2),
                  mean_degree=float(deg.mean()), degree_counts={int(k): int(v) for k, v in
                                                                zip(*np.unique(deg, return_counts=True))})
print("graph", graph_info)

labels = a.obs.domain.astype(str).values
cats = sorted(set(labels), key=int)
code = np.array([cats.index(l) for l in labels])
r, c = sp.triu(A, k=1).nonzero()


def same_frac(cd, rr=r, cc=c):
    return float(np.mean(cd[rr] == cd[cc]))


def per_domain_self(cd, rr=r, cc=c):
    out = np.zeros(len(cats))
    for k in range(len(cats)):
        m = (cd[rr] == k) | (cd[cc] == k)
        out[k] = np.mean((cd[rr] == k) & (cd[cc] == k) & m) / max(np.mean(m), 1e-12)
    return out


obs_global = same_frac(code)
obs_dom = per_domain_self(code)

# ---------------- Null A: label permutation on fixed graph ----------------
rng = np.random.default_rng(SEED)
nullA = np.empty(N_PERM)
nullA_dom = np.empty((N_PERM, len(cats)))
for i in range(N_PERM):
    p = rng.permutation(code)
    nullA[i] = same_frac(p)
    nullA_dom[i] = per_domain_self(p)

# ---------------- Null B: coordinate permutation + graph rebuild ----------------
rng_b = np.random.default_rng(SEED + 1)
# Uses the same coordinates (obsm["spatial"]) as the observed graph. Note: because spots are
# permuted among a fixed set of positions, the rebuilt graph is the observed graph with relabelled
# nodes, so Null B is equal in distribution to Null A. It is kept as an implementation check of the
# graph-building step, not as an independent null.
pos = np.asarray(a.obsm["spatial"], dtype=float)
nullB = np.empty(N_PERM)
nullB_edges = np.empty(N_PERM, dtype=int)
tmp = sc.AnnData(np.zeros((a.n_obs, 1), dtype=np.float32))
for i in range(N_PERM):
    perm = rng_b.permutation(a.n_obs)                 # spot i moves to position perm[i]
    tmp.obsm["spatial"] = pos[perm]
    sq.gr.spatial_neighbors(tmp, coord_type="grid", n_neighs=6, n_rings=1)
    B = tmp.obsp["spatial_connectivities"]
    rb, cb = sp.triu(((B + B.T) > 0), k=1).nonzero()
    nullB[i] = same_frac(code, rb, cb)
    nullB_edges[i] = len(rb)
assert (nullB_edges == len(r)).all(), "Null B graph must have the same edge count as the observed graph"


def emp_p(obs, null):           # one-sided, H1: observed is MORE coherent than null
    return (1 + np.sum(null >= obs)) / (len(null) + 1)


res = pd.DataFrame([
    dict(statistic="global same-domain edge fraction", null="A: label permutation, fixed graph",
         observed=obs_global, null_mean=nullA.mean(), null_sd=nullA.std(ddof=1),
         null_max=nullA.max(), z=(obs_global - nullA.mean()) / nullA.std(ddof=1),
         empirical_p=emp_p(obs_global, nullA), n_perm=N_PERM, seed=SEED),
    dict(statistic="global same-domain edge fraction", null="B: coordinate permutation, graph rebuilt",
         observed=obs_global, null_mean=nullB.mean(), null_sd=nullB.std(ddof=1),
         null_max=nullB.max(), z=(obs_global - nullB.mean()) / nullB.std(ddof=1),
         empirical_p=emp_p(obs_global, nullB), n_perm=N_PERM, seed=SEED + 1),
])
# expected value under random labelling (analytic) = sum_k n_k(n_k-1) / (N(N-1))
nk = np.bincount(code); N = len(code)
res["analytic_expectation"] = float(np.sum(nk * (nk - 1)) / (N * (N - 1)))
res.to_csv(OUT / "table_t5a_global_statistic_vs_null.csv", index=False)
print(res.T)

dom = pd.DataFrame(dict(domain=cats, n_spots=nk, observed_self_fraction=obs_dom,
                        null_mean=nullA_dom.mean(0), null_sd=nullA_dom.std(0, ddof=1),
                        empirical_p=[emp_p(obs_dom[k], nullA_dom[:, k]) for k in range(len(cats))]))
dom["fold_over_null"] = dom.observed_self_fraction / dom.null_mean
dom.round(4).to_csv(OUT / "table_t2a_per_domain_self_neighbour.csv", index=False)
print(dom.round(3).to_string())

# ---------------- S3 neighbourhood enrichment (squidpy) ----------------
sq.gr.nhood_enrichment(a, cluster_key="domain", n_perms=N_PERM, seed=SEED, show_progress_bar=False)
z = pd.DataFrame(a.uns["domain_nhood_enrichment"]["zscore"], index=cats, columns=cats)
z.round(2).to_csv(OUT / "table_t2b_nhood_enrichment_z.csv")
fig, ax = plt.subplots(figsize=(8, 6.8))
sns.heatmap(z, cmap="RdBu_r", center=0, vmin=-30, vmax=30, ax=ax, cbar_kws={"label": "z-score (vs 1000 label perms)"})
ax.set_title("Neighbourhood enrichment of expression domains\n(6-NN hex grid; diagonal = self-adjacency)")
ax.set_xlabel("domain"); ax.set_ylabel("domain")
fig.tight_layout(); fig.savefig(OUT / "fig_t2a_nhood_enrichment.png", dpi=130); plt.close(fig)

# strongest off-diagonal contacts (which domains border which)
off = z.where(~np.eye(len(cats), dtype=bool)).stack().sort_values(ascending=False)
off = off[[i < j for i, j in zip(off.index.get_level_values(0).astype(int), off.index.get_level_values(1).astype(int))]]
off.head(10).rename("z").to_csv(OUT / "table_t2c_top_domain_contacts.csv")
print("top contacts\n", off.head(8).round(1))

# ---------------- S4 Moran's I of domain indicators ----------------
ind = pd.get_dummies(a.obs.domain.astype(str)).astype(float)[cats]
ad_ind = sc.AnnData(ind.values, obs=a.obs[[]].copy())
ad_ind.var_names = [f"domain_{k}" for k in cats]
ad_ind.obsp["spatial_connectivities"] = A.astype(float)
sq.gr.spatial_autocorr(ad_ind, mode="moran", connectivity_key="spatial_connectivities",
                       n_perms=N_PERM, seed=SEED, show_progress_bar=False, n_jobs=1)
mi = ad_ind.uns["moranI"][["I", "pval_norm", "pval_sim", "pval_z_sim"]]
mi.round(4).to_csv(OUT / "table_t2d_domain_indicator_moransI.csv", index_label="domain")
print(mi.round(3))

# ---------------- T5 figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4))
for ax, null, ttl in [(axes[0], nullA, "Null A: permute labels on fixed graph"),
                      (axes[1], nullB, "Null B: permute positions, rebuild graph\n(implementation check; equivalent to A)")]:
    ax.hist(null, bins=40, color="grey", label=f"null (n={N_PERM})")
    ax.axvline(obs_global, color="crimson", lw=2, label=f"observed = {obs_global:.3f}")
    ax.set_xlabel("same-domain edge fraction"); ax.set_ylabel("count")
    ax.set_title(f"{ttl}\nnull {null.mean():.3f} ± {null.std(ddof=1):.3f}; p_emp = {emp_p(obs_global, null):.4f}")
    ax.legend(fontsize=8)
ax = axes[2]
o = np.argsort(-dom.fold_over_null.values)
ax.bar(range(len(cats)), dom.observed_self_fraction.values[o], color="crimson", label="observed")
ax.errorbar(range(len(cats)), dom.null_mean.values[o], yerr=2 * dom.null_sd.values[o], fmt="o",
            color="k", ms=3, label="null mean ± 2 SD")
ax.set_xticks(range(len(cats)), [cats[i] for i in o]); ax.set_xlabel("domain")
ax.set_ylabel("self-neighbour fraction"); ax.set_title("Per-domain self-adjacency vs null A"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "fig_t5_observed_vs_null.png", dpi=130); plt.close(fig)

json.dump(dict(graph="squidpy spatial_neighbors coord_type=grid n_neighs=6 n_rings=1 on obsm['spatial'] (Visium hex lattice)",
               graph_info=graph_info, n_perm=N_PERM, seed_nullA=SEED, seed_nullB=SEED + 1,
               nhood_enrichment_seed=SEED, tail="one-sided (observed >= null)"),
          open(OUT / "spatial_stats_settings.json", "w"), indent=2)
a.write_h5ad(OUT / "visium_domains.h5ad", compression="gzip")
