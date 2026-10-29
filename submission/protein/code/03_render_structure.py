"""P4: annotated EGFR structure (4HJO chain A, erlotinib) with all three
variant sites, labelled in CANONICAL numbering (author = canonical - 24).

Outputs: fig_p4a_4hjo_overview.png, fig_p4b_4hjo_atp_site.png (PyMOL ray-traced),
         fig_p4_annotated_structure.png (both panels + legend), egfr_4hjo_interactive.html (py3Dmol).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pymol import cmd

HERE = Path(__file__).resolve()
DATA = HERE.parents[3] / "data" / "protein"
OUT = HERE.parents[1] / "outputs"
OFF = 24  # 4HJO author numbering = canonical - 24 (verified in 01_map_and_measure.py)


def a(canon):
    return canon - OFF


REG = {  # canonical ranges -> colours
    "P-loop 718-726": ((718, 726), "orange"),
    "alphaC helix 752-769": ((752, 769), "yellow"),
    "hinge 790-791": ((790, 791), "cyan"),
    "activation segment 855-884": ((855, 884), "violet"),
}
SITES = {"G719": 719, "T790": 790, "L858": 858}

cmd.reinitialize()
cmd.load(str(DATA / "4hjo.pdb"), "egfr")
cmd.remove("solvent or hydro")
cmd.hide("everything")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.show("cartoon", "egfr and polymer")
cmd.color("grey80", "egfr and polymer")
for name, ((lo, hi), col) in REG.items():
    cmd.color(col, f"egfr and polymer and resi {a(lo)}-{a(hi)}")
cmd.show("sticks", "resn AQ4"); cmd.color("green", "resn AQ4 and elem C"); cmd.util.cnc("resn AQ4")
for lab, p in SITES.items():
    sel = f"site_{lab}"
    cmd.select(sel, f"egfr and polymer and resi {a(p)}")
    cmd.show("spheres", f"{sel} and not name N+C+O")
    cmd.set("sphere_scale", 0.6, sel)
    cmd.color("red", f"{sel} and elem C")
    cmd.label(f"{sel} and name CA", f'"{lab}"')
    cmd.set("label_position", {"G719": (4.0, -3.0, 25), "T790": (3.5, -2.5, 25), "L858": (-5.0, -2.5, 25)}[lab],
            f"{sel} and name CA")
cmd.select("v948r", f"egfr and resi {a(948)}")
cmd.show("sticks", "v948r and not name N+C+O"); cmd.color("blue", "v948r and elem C")
cmd.label("v948r and name CA", '"V948R"'); cmd.set("label_position", (0.0, 4.0, 25), "v948r and name CA")
cmd.select("cat", f"egfr and resi {a(745)}+{a(762)}+{a(837)}+{a(855)} and not name N+C+O")
cmd.show("sticks", "cat"); cmd.color("slate", "cat and elem C")
cmd.set("label_size", 22); cmd.set("label_color", "black"); cmd.set("label_font_id", 7)
cmd.set("label_outline_color", "white"); cmd.set("float_labels", 1)
cmd.set("cartoon_transparency", 0.15)
cmd.set("antialias", 2); cmd.set("ray_shadows", 0); cmd.set("depth_cue", 0)
cmd.orient("egfr and polymer")
cmd.turn("x", -60); cmd.turn("y", 20)
cmd.zoom("egfr and polymer", -6)
cmd.png(str(OUT / "fig_p4a_4hjo_overview.png"), width=1400, height=1100, dpi=200, ray=1)

# close-up of ATP site with measured distances (min heavy-atom, from 01 script)
cmd.label("all", "''")
cmd.hide("sticks", "cat or v948r")
cmd.show("sticks", f"egfr and resi {a(745)}+{a(762)} and not name N+C+O")
offsets = {"G719": (2.5, -1.5, 3), "T790": (2.0, 1.5, 3), "L858": (-3.5, -1.0, 3)}
for lab, p in SITES.items():
    cmd.label(f"site_{lab} and name CA", f'"{lab}"')
    cmd.set("label_position", offsets[lab], f"site_{lab} and name CA")
cmd.label(f"egfr and resi {a(745)} and name NZ", '"K745"'); cmd.set("label_position", (3.5, 1.0, 3), f"egfr and resi {a(745)} and name NZ")
cmd.label(f"egfr and resi {a(762)} and name CD", '"E762"'); cmd.set("label_position", (-1.5, 1.5, 3), f"egfr and resi {a(762)} and name CD")
cmd.distance("d790", "site_T790 and not name N+C+O+CA", "resn AQ4", cutoff=3.6, mode=0)
cmd.distance("d719", "site_G719", "resn AQ4", cutoff=5.1, mode=0)
cmd.distance("dKE", f"egfr and resi {a(745)} and name NZ", f"egfr and resi {a(762)} and name OE1+OE2", mode=0, cutoff=14)
for dname, col in [("d790", "red"), ("d719", "red"), ("dKE", "blue")]:
    cmd.color(col, dname)
cmd.hide("labels", "d*")
cmd.set("dash_width", 4); cmd.set("label_size", 30)
cmd.set("cartoon_transparency", 0.55, "egfr and polymer and not resi " + f"{a(718)}-{a(726)}+{a(752)}-{a(769)}+{a(855)}-{a(884)}")
focus = f"resn AQ4 or site_G719 or site_T790 or site_L858 or (egfr and resi {a(745)}+{a(762)})"
cmd.orient(focus)
cmd.zoom(focus, 2)
cmd.png(str(OUT / "fig_p4b_4hjo_atp_site.png"), width=1400, height=1100, dpi=200, ray=1)

# compose
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
for ax, f, t in zip(axes, ["fig_p4a_4hjo_overview.png", "fig_p4b_4hjo_atp_site.png"],
                    ["A. EGFR kinase domain, PDB 4HJO chain A\n(inactive, erlotinib)",
                     "B. ATP site close-up (red dashes = closest site–erlotinib distance;\nblue = broken K745–E762 salt bridge)"]):
    ax.imshow(plt.imread(OUT / f)); ax.axis("off"); ax.set_title(t, fontsize=12)
handles = [mpatches.Patch(color=c, label=l) for l, c in [
    ("P-loop 718-726 (UniProt ATP-binding)", "#ffa500"), ("αC helix 752-769", "#ffff00"),
    ("hinge 790-791", "#00cccc"), ("activation segment 855-884", "#ee82ee"),
    ("variant sites G719 / T790 / L858 (red spheres)", "red"), ("erlotinib (AQ4, green sticks)", "green"),
    ("K745, E762, D837, D855 (slate sticks)", "#7080ff"), ("V948R engineered (blue)", "blue")]]
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.025))
fig.text(0.5, 0.005, "Numbering: canonical UniProt P00533 (4HJO author number + 24). Distances (Table P2a): T790-erlotinib 3.48 Å, "
         "G719-erlotinib 5.06 Å, L858-D855 2.96 Å, K745 NZ-E762 OE 13.3 Å (salt bridge broken = αC-out, inactive). "
         "Minimum heavy-atom distances.", ha="center", fontsize=8.5)
fig.tight_layout(rect=(0, 0.10, 1, 0.95))
fig.savefig(OUT / "fig_p4_annotated_structure.png", dpi=160)

# interactive HTML (py3Dmol)
import py3Dmol
v = py3Dmol.view(width=900, height=650)
v.addModel(open(DATA / "4hjo.pdb").read(), "pdb")
v.setStyle({"cartoon": {"color": "lightgrey"}})
for name, ((lo, hi), col) in REG.items():
    v.setStyle({"resi": f"{a(lo)}-{a(hi)}"}, {"cartoon": {"color": {"yelloworange": "gold"}.get(col, col)}})
v.addStyle({"resn": "AQ4"}, {"stick": {"colorscheme": "greenCarbon"}})
for lab, p in SITES.items():
    v.addStyle({"resi": str(a(p))}, {"sphere": {"color": "red", "scale": 0.7}})
    v.addResLabels({"resi": str(a(p))}, {"fontSize": 14, "showBackground": True})
    v.addLabel(f"{lab} (canonical; 4HJO {a(p)})", {"fontSize": 12}, {"resi": str(a(p))})
v.zoomTo({"resn": "AQ4"})
html = v._make_html()
(OUT / "egfr_4hjo_interactive.html").write_text(html)
print("written")
