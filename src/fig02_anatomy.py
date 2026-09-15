"""Figure 2. The anatomy of a definition: what each one actually counts.

Every alternative definition is one row. Left block: one column per ICD-10 code used anywhere in
that condition's definitions, in code order; a filled tile means the definition counts that code
in hospital or death records, and a tile with a white notch cut into its corner also counts it
from primary care. Right block: the rules codes cannot show,
  drug   counts a drug purchase class (ATC A06A, laxatives)
  reimb  counts a drug reimbursement right
  excl   case must not also meet another definition (e.g. focal and not generalized)
  mode   the code must be the person's most frequent diagnosis of that kind
  ctrl   the control rule differs from the set's first definition (which people are excluded from
         controls; it can be more or fewer exclusions)
Two definitions with identical code rows can still differ in the right block, which is exactly
where some of them differ. Tiles show ICD-10 only; ICD-9 and ICD-8 rules are in the endpoint file.

The idea follows Dear Data's habit of layering several attributes onto a small hand-built mark
with one shared key, instead of a Venn diagram of code sets.

Sources: results/definitions/code_membership.tsv and endpoint_anatomy.tsv (src/02).
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import pandas as pd

import viz_style as vs
from itoju_labels import ENDPOINT_LABEL

ROOT = vs.ROOT
SETS = [
    ("Epilepsy", ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"]),
    ("Sleep apnoea", ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"]),
    ("Insomnia", ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"]),
    ("Constipation", ["K11_CONSTIPATION", "K11_OTHFUNC"]),
    ("ADHD", ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"]),
]
RULES = ["drug", "reimb", "excl", "mode", "ctrl"]


def main():
    vs.apply()
    mem = pd.read_csv(ROOT / "results" / "definitions" / "code_membership.tsv", sep="\t")
    ana = pd.read_csv(ROOT / "results" / "definitions" / "endpoint_anatomy.tsv", sep="\t").set_index("endpoint")

    n_rows = sum(len(d) for _, d in SETS) + len(SETS)
    fig = plt.figure(figsize=(vs.DOUBLE, 0.19 * n_rows + 0.9))
    ax = fig.add_axes([0.20, 0.09, 0.79, 0.84])
    ink, soft = vs.INK_2, vs.GRID
    y = 0
    max_x = 0
    for name, defs in SETS:
        codes = sorted(set(mem[mem.endpoint.isin(defs)].icd10))
        # collapse a three-character parent into its children when both are present, to keep tiles readable
        codes = [c for c in codes if not (len(c) == 3 and any(o.startswith(c) and len(o) > 3 for o in codes))] or codes
        ax.text(-0.6, y, name, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK)
        for k, c in enumerate(codes):
            ax.text(k, y, f"{c[:3]}.{c[3:]}" if len(c) > 3 else c, rotation=90, ha="center", va="center",
                    fontsize=7, color=vs.MUTED)
        rx0 = len(codes) + 1.2
        for j, rule in enumerate(RULES):
            ax.text(rx0 + j * 1.25, y, rule, ha="center", va="center", fontsize=7, color=vs.MUTED, rotation=90)
        y += 1.3
        ref_ctrl = str(ana.loc[defs[0], "control_exclude"]) + str(ana.loc[defs[0], "control_conditions"])
        for e in defs:
            a = ana.loc[e]
            ax.text(-0.6, y, ENDPOINT_LABEL[e], ha="right", va="center", fontsize=7, color=vs.INK_2)
            m = mem[mem.endpoint == e]
            hd = set(m[m.source.isin(["HD", "COD"])].icd10)
            op = set(m[m.source == "OUTPAT"].icd10)
            for k, c in enumerate(codes):
                covered = c in hd or any(h.startswith(c) for h in hd)
                ax.add_patch(Rectangle((k - 0.42, y - 0.36), 0.84, 0.72, facecolor=ink if covered else "white",
                                       edgecolor=ink if covered else soft, linewidth=0.5))
                if c in op or any(o.startswith(c) for o in op):
                    # a white notch cut into the tile: no trait color is borrowed for a non-trait meaning
                    ax.add_patch(Polygon([(k + 0.42, y - 0.36), (k + 0.42, y + 0.02), (k + 0.04, y - 0.36)],
                                         closed=True, facecolor=vs.SURFACE, edgecolor="none", zorder=3))
            sources = str(a.sources).split()
            flags = {
                "drug": "KELA_ATC" in sources,
                "reimb": "KELA_REIMB" in sources,
                "excl": isinstance(a.case_conditions, str) and a.case_conditions.strip() != "",
                "mode": bool(a.mode_rule),
                "ctrl": (str(a.control_exclude) + str(a.control_conditions)) != ref_ctrl,
            }
            for j, rule in enumerate(RULES):
                x = rx0 + j * 1.25
                if flags[rule]:
                    ax.scatter(x, y, s=34, marker="o", facecolor=ink, edgecolor=ink, zorder=3)
                else:
                    ax.scatter(x, y, s=10, marker="o", facecolor="white", edgecolor=soft, linewidths=0.6, zorder=3)
            y += 1
            max_x = max(max_x, rx0 + len(RULES) * 1.25)
        y += 0.6
    ax.set_xlim(-0.6, max_x)
    ax.set_ylim(y - 0.4, -0.8)
    ax.set_aspect("equal")
    ax.axis("off")
    # key, two rows
    kx, ky, ky2 = 0.20, 0.045, 0.015
    fig.add_artist(Rectangle((kx, ky - 0.008), 0.012, 0.018, transform=fig.transFigure, facecolor=ink, edgecolor=ink))
    fig.text(kx + 0.018, ky, "code counted (hospital or death records)", fontsize=7, color=vs.INK_2, va="center")
    fig.add_artist(Rectangle((kx + 0.33, ky - 0.008), 0.012, 0.018, transform=fig.transFigure, facecolor=ink,
                             edgecolor=ink))
    fig.add_artist(Polygon([(kx + 0.342, ky - 0.008), (kx + 0.342, ky + 0.003), (kx + 0.334, ky - 0.008)], closed=True,
                           transform=fig.transFigure, facecolor=vs.SURFACE, edgecolor="none"))
    fig.text(kx + 0.348, ky, "notched: also counted from primary care", fontsize=7, color=vs.INK_2, va="center")
    fig.add_artist(plt.Line2D([kx + 0.006], [ky2], marker="o", markersize=5, color=ink, ls="", transform=fig.transFigure))
    fig.text(kx + 0.018, ky2, "rule applies: drug purchase, reimbursement, exclusion of another definition, mode rule, "
             "control rule differs from the first row", fontsize=7, color=vs.INK_2, va="center")
    vs.save(fig, "fig02_anatomy")


if __name__ == "__main__":
    main()
