"""P1 + P2: verify EGFR variants, map canonical positions to every structure
by sequence alignment, and compute reproducible structural measurements.

Measurement rules (all heavy atoms; hydrogens/waters ignored; Angstrom):
  * min_dist(X, Y) = smallest Euclidean distance between any heavy atom of X
    and any heavy atom of Y.
  * Ligand = erlotinib (AQ4) in 4HJO/1M17, AMP-PNP (ANP) in 2ITX/2ITN/2ITV.
  * Neighbours = residues with any heavy atom <= 5.0 A of the site residue.
  * Contacts = number of heavy-atom pairs <= 4.5 A.
  * Relative SASA = freesasa (Lee-Richards, probe 1.4 A) on the protein chain
    alone (ligand and waters removed), relativeTotal.
  * Relative B = mean B of residue / mean B of all protein heavy atoms in chain.
"""
from structlib import *  # noqa: F401,F403  (constants + helpers)
import freesasa

# ------------------------------------------------------------------ P1 ----
checks = []
for name, (pos, wt, mt) in VARIANTS.items():
    checks.append(dict(variant=name, canonical_pos=pos, stated_wt=wt, uniprot_residue=SEQ[pos - 1],
                       wt_matches=SEQ[pos - 1] == wt, context=SEQ[pos - 6:pos - 1] + "[" + SEQ[pos - 1] + "]" + SEQ[pos:pos + 5]))
checks = pd.DataFrame(checks)
print(checks.to_string())
assert checks.wt_matches.all()
motifs = {m: SEQ.find(m, 700) + 1 for m in ["GSGAFG", "HRDLAARN", "DFG", "ALE"]}  # EGFR has ALE in place of APE
print("motifs (1-based canonical start):", motifs)

mapping_rows, meas_rows, neigh_rows, state_rows = [], [], [], []
for sid, (ch, lig, desc) in STRUCTS.items():
    st = load(sid)
    chain = st[0][ch]
    m, mism, res = map_chain(st, chain)
    inv = {v: k for k, v in m.items()}
    offs = sorted({v - k for k, v in m.items()})
    cov = (min(m.values()), max(m.values()))
    prot = [r for r in res]
    ligres = [r for r in chain if r.name == lig] if lig else []
    if lig and not ligres:   # ligand may sit in its own chain
        ligres = [r for c in st[0] for r in c if r.name == lig][:1]
    allB = np.mean([a.b_iso for r in prot for a in r])

    # SASA on protein-only chain
    tmp = gemmi.Structure(); tmp.add_model(gemmi.Model("1")); nc = gemmi.Chain(ch)
    for r in prot:
        nc.add_residue(r)
    tmp[0].add_chain(nc); tmp.setup_entities()
    tmp_pdb = OUT / f"_tmp_{sid}.pdb"; tmp.write_pdb(str(tmp_pdb))
    fs = freesasa.Structure(str(tmp_pdb)); sres = freesasa.calc(fs).residueAreas()[ch]
    tmp_pdb.unlink()

    def R(canon):
        a = inv.get(canon)
        return None if a is None else next(r for r in prot if r.seqid.num == a)

    alphaC = [R(i) for i in range(752, 770) if R(i) is not None]
    # conformation markers
    k, e = R(745), R(762)
    ke = mind(heavy([k]) if k else np.empty((0, 3)), heavy([e]) if e else np.empty((0, 3)))
    k_nz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in k if a.name == "NZ"]) if k else np.empty((0, 3))
    e_cd = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in e if a.name in ("OE1", "OE2")]) if e else np.empty((0, 3))
    state_rows.append(dict(structure=sid.upper(), description=desc,
                           K745NZ_E762OE_min=mind(k_nz, e_cd),
                           conformation=("alphaC-in (K745-E762 salt bridge)" if mind(k_nz, e_cd) < 4.0
                                         else "alphaC-out (salt bridge broken)")))
    missing_sites = [p for p in (719, 790, 858) if p not in inv]
    mapping_rows.append(dict(structure=sid.upper(), chain=ch, description=desc,
                             author_minus_canonical=",".join(str(-o) for o in offs),
                             resolved_canonical_range=f"{cov[0]}-{cov[1]}", n_resolved=len(m),
                             **{f"{v}_author_res": (f"{R(p).name}{R(p).seqid.num}" if R(p) else "UNRESOLVED")
                                for v, (p, _, _) in VARIANTS.items()},
                             sequence_mismatches_vs_UniProt="; ".join(mism) or "none",
                             missing_sites=",".join(map(str, missing_sites)) or "none"))

    for v, (p, wt, mt) in VARIANTS.items():
        r = R(p)
        if r is None:
            meas_rows.append(dict(structure=sid.upper(), variant=v, resolved=False)); continue
        X = heavy([r])
        row = dict(structure=sid.upper(), variant=v, resolved=True, residue_in_structure=f"{r.name}{r.seqid.num}",
                   min_dist_to_ligand=mind(X, heavy(ligres)) if ligres else np.nan,
                   ligand=lig or "none",
                   min_dist_to_alphaC=mind(X, heavy([a for a in alphaC if a is not r])),
                   contacts_alphaC_4p5A=int(sum((np.linalg.norm(X[:, None] - heavy([a])[None], axis=2) <= 4.5).sum()
                                                 for a in alphaC if a is not r)),
                   rel_SASA=sres[str(r.seqid.num)].relativeTotal if str(r.seqid.num) in sres else np.nan,
                   rel_B=float(np.mean([a.b_iso for a in r]) / allB))
        for lab, (kp, names) in KEY_ATOMS.items():
            kr = R(kp)
            A = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in kr if a.name in names]) if kr else np.empty((0, 3))
            row[f"min_dist_to_{lab.split()[0]}"] = mind(X, A)
        # backbone torsions for the Gly site
        if p == 719:
            ph = gemmi.calculate_phi_psi(chain.previous_residue(r), r, chain.next_residue(r))
            row["phi"], row["psi"] = np.degrees(ph[0]), np.degrees(ph[1])
        meas_rows.append(row)
        # neighbours within 5 A (canonical numbering)
        for q in prot:
            if q is r:
                continue
            d = mind(X, heavy([q]))
            if d <= 5.0:
                neigh_rows.append(dict(structure=sid.upper(), variant=v, neighbour=f"{q.name}{m[q.seqid.num]}",
                                       min_dist=round(d, 2)))
        for L in ligres:
            d = mind(X, heavy([L]))
            if d <= 5.0:
                neigh_rows.append(dict(structure=sid.upper(), variant=v, neighbour=f"ligand {L.name}", min_dist=round(d, 2)))

mp = pd.DataFrame(mapping_rows); ms = pd.DataFrame(meas_rows); nb = pd.DataFrame(neigh_rows); cs = pd.DataFrame(state_rows)
checks.to_csv(OUT / "table_p1a_variant_sequence_check.csv", index=False)
mp.to_csv(OUT / "table_p1b_structure_mapping.csv", index=False)
ms.round(2).to_csv(OUT / "table_p2a_site_measurements.csv", index=False)
nb.to_csv(OUT / "table_p2b_neighbours_within_5A.csv", index=False)
cs.round(2).to_csv(OUT / "table_p2c_conformation_state.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(mp.drop(columns=["description"]).to_string())
print(cs.round(2).to_string())
print(ms.round(2).to_string())

# ------------------------------------------------- pLDDT (AlphaFold) ----
conf = json.load(open(DATA / "AF-P00533-F1-confidence_v6.json"))
pl = dict(zip(conf["residueNumber"], conf["confidenceScore"]))
plddt = pd.DataFrame([dict(variant=v, canonical_pos=p, pLDDT=pl[p],
                           window_mean_pLDDT_pm5=np.mean([pl[i] for i in range(p - 5, p + 6)])) for v, (p, _, _) in VARIANTS.items()])
plddt["kinase_domain_mean_pLDDT"] = np.mean([pl[i] for i in range(712, 980)])
plddt.round(1).to_csv(OUT / "table_p1c_alphafold_plddt.csv", index=False)
print(plddt.round(1))

# ------------------------------------------- paralog conservation ----
para = {}
for rec in open(DATA / "erbb_paralogs.fasta").read().split(">")[1:]:
    head, *seq = rec.splitlines(); para[head.split("|")[2].split()[0]] = "".join(seq)
cons = []
for nm, s in para.items():
    aln = aligner.align(SEQ, s)[0]
    mpos = {i + 1: j + 1 for (a0, a1), (b0, b1) in zip(*aln.aligned) for i, j in zip(range(a0, a1), range(b0, b1))}
    for v, (p, _, _) in VARIANTS.items():
        j = mpos.get(p)
        cons.append(dict(paralog=nm, variant=v, paralog_pos=j, paralog_residue=s[j - 1] if j else "-",
                         local=s[j - 4:j + 3] if j else ""))
cons = pd.DataFrame(cons)
cons.to_csv(OUT / "table_p3a_paralog_conservation.csv", index=False)
print(cons.pivot(index="variant", columns="paralog", values="paralog_residue"))

# AlphaMissense (retrieved, not computed)
am = pd.read_csv(DATA / "AF-P00533-F1-aa-substitutions.csv")
am = am[am.protein_variant.isin(list(VARIANTS))]
am.to_csv(OUT / "table_p3b_alphamissense.csv", index=False)
print(am)

# ------------------------------- engineered V948R distance (4HJO) ----
st = load("4hjo"); ch = st[0]["A"]
m, _, res = map_chain(st, ch)
byc = {m[r.seqid.num]: r for r in res}
aq4 = [r for r in ch if r.name == "AQ4"]
v948 = pd.DataFrame([dict(target=t, min_dist_from_R948=round(mind(heavy([byc[948]]), heavy(x)), 1))
                     for t, x in [("G719", [byc[719]]), ("T790", [byc[790]]), ("L858", [byc[858]]), ("erlotinib AQ4", aq4)]])
v948.to_csv(OUT / "table_p2f_V948R_distances_4hjo.csv", index=False)
print(v948)
