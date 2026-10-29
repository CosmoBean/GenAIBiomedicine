"""P2/P3: variant-specific structural checks.

(a) T790M  - superpose the experimental Met790 structures (2JIT, 4I22) onto
             the erlotinib complexes (4HJO inactive, 1M17 active) and the WT
             AMP-PNP complex (2ITX) by kinase-domain CA atoms (canonical
             712-979, identical-residue pairs only), then measure Met790 vs
             Thr790 side-chain distance to the ligand.
(b) L858R  - environment of residue 858 in inactive (4HJO) vs active-like
             (1M17, 2ITX) WT structures: buried hydrophobic neighbours and
             relative SASA; and Arg858 environment in 2ITV.
(c) G719S  - Ser719 OG to AMP-PNP phosphate oxygens in 2ITN, vs the Gly719
             CA in 2ITX; P-loop CA RMSD 2ITN vs 2ITX after superposition.
"""
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
OUT = HERE.parents[1] / "outputs"
DATA = HERE.parents[3] / "data" / "protein"

from structlib import load, map_chain, heavy, mind  # noqa: E402

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR", "CYS"}


def chain_info(sid, ch="A"):
    st = load(sid)
    c = st[0][ch]
    m, _, res = map_chain(st, c)
    by_can = {m[r.seqid.num]: r for r in res}
    return st, c, by_can


def ca_dict(by_can, lo=712, hi=979):
    return {k: r["CA"][0].pos for k, r in by_can.items() if lo <= k <= hi and r.find_atom("CA", "*")}


def superpose(mov_bycan, ref_bycan, positions=None):
    a, b = ca_dict(mov_bycan), ca_dict(ref_bycan)
    keys = sorted(set(a) & set(b))
    if positions is not None:
        keys = [k for k in keys if k in positions]
    sup = gemmi.superpose_positions([b[k] for k in keys], [a[k] for k in keys])
    return sup, len(keys)


def coords(res, names=None):
    return np.array([[at.pos.x, at.pos.y, at.pos.z] for at in res if names is None or at.name in names])


def ligand(st, name):
    return [r for c in st[0] for r in c if r.name == name]


rows = []
# ---------------- (a) T790M ----------------
refs = {"4hjo": "AQ4", "1m17": "AQ4", "2itx": "ANP"}
for mut in ["2jit", "4i22"]:
    _, _, mb = chain_info(mut)
    met = mb[790]
    for ref, lig in refs.items():
        rst, _, rb = chain_info(ref)
        L = heavy(ligand(rst, lig)[:1])
        thr = coords(rb[790], ["OG1", "CG2", "CB"])
        # pocket = reference residues with any heavy atom within 10 A of the ligand
        pocket = {k for k, r in rb.items() if mind(heavy([r]), L) <= 10.0}
        for mode, pos in [("global kinase-domain CA", None), ("local ATP-pocket CA (<=10 A of ligand)", pocket)]:
            sup, n = superpose(mb, rb, pos)
            sc_atoms = np.array([list(sup.transform.apply(at.pos).tolist()) for at in met if at.name in ("CG", "SD", "CE")])
            rows.append(dict(check="T790M side chain vs ligand (superposed)", superposition=mode,
                             mutant_structure=mut.upper(), reference=f"{ref.upper()} ({lig})",
                             n_CA_superposed=n, CA_RMSD=round(sup.rmsd, 2),
                             WT_Thr790_sidechain_to_ligand=round(mind(thr, L), 2),
                             Met790_sidechain_to_ligand=round(mind(sc_atoms, L), 2),
                             n_ligand_atoms_within_3A_of_Met=int((np.linalg.norm(L[:, None] - sc_atoms[None], axis=2) < 3.0).any(1).sum())))

# ---------------- (c) G719S ----------------
st_n, _, nb = chain_info("2itn")
st_x, _, xb = chain_info("2itx")
anp_n = ligand(st_n, "ANP")[0]
phos_O = coords(anp_n, [a.name for a in anp_n if a.name.startswith("O") and a.name[1] in "123" and a.name[-1] in "ABG"])
og = coords(nb[719], ["OG"])
sup, n = superpose(nb, xb)
ploop = [k for k in range(718, 727) if k in nb and k in xb]
d = [np.linalg.norm(np.array(sup.transform.apply(nb[k]["CA"][0].pos).tolist()) - np.array(xb[k]["CA"][0].pos.tolist()))
     for k in ploop]
anp_x = ligand(st_x, "ANP")[0]
phos_Ox = coords(anp_x, [a.name for a in anp_x if a.name.startswith("O") and a.name[1] in "123" and a.name[-1] in "ABG"])
rows.append(dict(check="G719S Ser OG to AMP-PNP phosphate O", superposition="global kinase-domain CA", mutant_structure="2ITN", reference="2ITX (WT)",
                 n_CA_superposed=n, CA_RMSD=round(sup.rmsd, 2),
                 WT_Gly719_CA_to_phosphateO=round(mind(coords(xb[719], ["CA"]), phos_Ox), 2),
                 Ser719_OG_to_phosphateO=round(mind(og, phos_O), 2),
                 Ploop_718_726_CA_RMSD=round(float(np.sqrt(np.mean(np.square(d)))), 2),
                 Ploop_max_CA_shift=round(float(np.max(d)), 2)))

# ---------------- (b) L858 environment ----------------
env = []
for sid in ["4hjo", "1m17", "2itx", "2itv", "4i22"]:
    st, c, b = chain_info(sid)
    r = b[858]
    X = heavy([r])
    nbrs = []
    for k, q in b.items():
        if q is r or abs(k - 858) <= 1:
            continue
        if mind(X, heavy([q])) <= 4.5:
            nbrs.append((k, q.name))
    hyd = [f"{n[:1]}{n[1:].lower()}{k}" for k, n in nbrs if n in HYDROPHOBIC]
    acid = [f"{n}{k}" for k, n in nbrs if n in ("ASP", "GLU")]
    env.append(dict(structure=sid.upper(), residue858=r.name, n_neighbours_4p5A=len(nbrs),
                    hydrophobic_neighbours=", ".join(hyd), n_hydrophobic=len(hyd),
                    acidic_neighbours=", ".join(acid),
                    alphaC_neighbours=", ".join(f"{n}{k}" for k, n in nbrs if 752 <= k <= 769)))
env = pd.DataFrame(env)
env.to_csv(OUT / "table_p2e_L858_environment.csv", index=False)

res = pd.DataFrame(rows)
res.to_csv(OUT / "table_p2d_T790M_G719S_checks.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(res.drop(columns=[c for c in res if c.startswith(("WT_Gly", "Ser719", "Ploop"))]).iloc[:-1].to_string()); print(res.iloc[-1].dropna().to_string()); print(env.to_string())
