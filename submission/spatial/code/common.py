"""Shared settings for the Visium mouse-brain analysis (Task 3)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "visium"
COUNTS = "V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5"
OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(exist_ok=True)

SEED = 0
QC = dict(min_counts=5000, min_genes=1000, max_pct_mt=30.0, min_cells=10)
N_HVG = 2000          # for PCA / clustering
N_PCS = 30
N_NEIGHBORS = 15
RES = 0.5
N_SVG_GENES = 200     # gene set tested for spatial autocorrelation + DE (T3/T4)
N_PERM = 1000         # randomizations for T5 (>=100 required)
