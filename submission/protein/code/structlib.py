"""Shared EGFR data, constants and structure helpers for protein/ scripts."""

import json
from pathlib import Path

import freesasa
import gemmi
import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve()
DATA = HERE.parents[3] / "data" / "protein"
OUT = HERE.parents[1] / "outputs"
OUT.mkdir(exist_ok=True)

uni = json.load(open(DATA / "P00533.json"))
SEQ = uni["sequence"]["value"]
VARIANTS = {"G719S": (719, "G", "S"), "T790M": (790, "T", "M"), "L858R": (858, "L", "R")}

# Regions (canonical numbering). Source given for each.
REGIONS = {
    "Protein kinase domain (UniProt)": (712, 979),
    "Gly-rich P-loop / ATP-binding (UniProt binding site)": (718, 726),
    "alphaC helix (HELIX records 4HJO 755-769, 2ITX 752-768)": (752, 769),
    "Hinge / ATP-binding (UniProt 790-791)": (790, 791),
    "Activation segment DFG...ALE (sequence motifs)": (855, 884),
    "A-loop short helix in inactive 4HJO (HELIX record)": (857, 862),
}
KEY_ATOMS = {"K745 (catalytic Lys)": (745, ["NZ"]), "E762 (alphaC Glu)": (762, ["CD", "OE1", "OE2"]),
             "D837 (HRD catalytic Asp)": (837, ["CG", "OD1", "OD2"]),
             "D855 (DFG Asp)": (855, ["CG", "OD1", "OD2"])}

STRUCTS = {  # id: (chain, ligand, description)
    "4hjo": ("A", "AQ4", "WT (+V948R) inactive, erlotinib, 2.75 A [primary]"),
    "1m17": ("A", "AQ4", "WT active-like, erlotinib, 2.6 A"),
    "2itx": ("A", "ANP", "WT, AMP-PNP, 2.98 A"),
    "2itn": ("A", "ANP", "G719S, AMP-PNP, 2.47 A"),
    "2itv": ("A", "ANP", "L858R, AMP-PNP, 2.47 A"),
    "2jit": ("A", None, "T790M apo, 3.1 A"),
    "4i22": ("A", "IRE", "L858R+T790M (+V948R), gefitinib, 1.71 A"),
    "AF-P00533-F1-model_v6": ("A", None, "AlphaFold DB v6 model (full length)"),
}

aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")


def load(sid):
    st = gemmi.read_structure(str(DATA / f"{sid}.pdb"))
    st.remove_hydrogens(); st.remove_waters()
    return st


def map_chain(st, chain):
    """Map author residue numbers to canonical UniProt numbers.

    The full deposited construct sequence (SEQRES, incl. unresolved residues)
    is aligned to UniProt, and each resolved residue is placed through its
    label_seq index. Aligning SEQRES instead of only the resolved residues
    keeps gap placement unambiguous around disordered loops.
    """
    st.setup_entities()
    st.assign_label_seq_id(force=True)
    ent = st.get_entity_of(chain.get_polymer())
    full = gemmi.one_letter_code(ent.full_sequence).upper()
    aln = aligner.align(full, SEQ)[0]
    lab2can = {i + 1: j + 1 for (a0, a1), (b0, b1) in zip(*aln.aligned)
               for i, j in zip(range(a0, a1), range(b0, b1))}
    res = [r for r in chain.get_polymer()]
    m, mism = {}, []
    for r in res:
        j = lab2can.get(r.label_seq)
        if j is None:
            continue
        m[r.seqid.num] = j
        aa = gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
        if aa != SEQ[j - 1]:
            mism.append(f"{r.name}{r.seqid.num}->canon {SEQ[j-1]}{j}")
    return m, mism, [r for r in res if r.seqid.num in m]


def heavy(res_list):
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for r in res_list for a in r if a.element.name != "H"])


def mind(A, B):
    if len(A) == 0 or len(B) == 0:
        return np.nan
    return float(np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)))


