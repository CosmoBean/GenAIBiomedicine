"""Extra control for T5: are the domains driven by spatially structured
technical covariates (sequencing depth, mitochondrial fraction)?

Re-run the identical expression-only pipeline but regress out log total
counts and % mito before scaling/PCA. Compare (i) domains with ARI and
(ii) the same global same-domain edge fraction + label-permutation null.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score as ARI

from common import OUT, SEED, N_PCS, N_NEIGHBORS, RES, N_PERM

a = sc.read_h5ad(OUT / "visium_domains.h5ad")
x = sc.AnnData(a.X.copy(), obs=a.obs[["total_counts", "pct_counts_mt"]].copy(), var=a.var[["highly_variable"]].copy())
x.obs["log_total_counts"] = np.log1p(x.obs.total_counts)
x = x[:, x.var.highly_variable].copy()
sc.pp.regress_out(x, ["log_total_counts", "pct_counts_mt"])
sc.pp.scale(x, max_value=10)
sc.tl.pca(x, n_comps=50, random_state=SEED)
sc.pp.neighbors(x, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)
sc.tl.leiden(x, resolution=RES, key_added="domain_reg", random_state=SEED, flavor="igraph",
             n_iterations=2, directed=False)

A = a.obsp["spatial_connectivities"]
r, c = sp.triu(((A + A.T) > 0), k=1).nonzero()
code = x.obs.domain_reg.cat.codes.to_numpy()
obs = float(np.mean(code[r] == code[c]))
rng = np.random.default_rng(SEED)
null = np.array([np.mean((p := rng.permutation(code))[r] == p[c]) for _ in range(N_PERM)])
out = pd.Series(dict(
    n_domains_original=a.obs.domain.nunique(), n_domains_depth_regressed=x.obs.domain_reg.nunique(),
    ARI_original_vs_regressed=ARI(a.obs.domain, x.obs.domain_reg),
    same_domain_edge_frac_regressed=obs, null_mean=null.mean(), null_sd=null.std(ddof=1),
    empirical_p=(1 + np.sum(null >= obs)) / (N_PERM + 1)))
out.round(4).to_csv(OUT / "table_t5c_depth_regressed_domains.csv", header=["value"])
print(out.round(4))
ct = pd.crosstab(a.obs.domain, x.obs.domain_reg)
ct.to_csv(OUT / "table_t5d_original_vs_depth_regressed_crosstab.csv")
