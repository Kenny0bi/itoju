"""Figure 2. The anatomy of a definition: what each one actually counts.

Every alternative definition is one row. Left block: one tile per ICD-10 code used anywhere in that
condition's definitions, in code order. An indigo tile counts the code in hospital or death records;
a teal corner means the code also counts from primary care. Right block, aligned for every condition:
the rules codes cannot show, each drawn as its own small glyph, in the manner of Dear Data's hand-built
marks with one shared key (right):
  pill        counts a drug purchase class (ATC A06A, laxatives)
  card        counts a drug reimbursement right
  slash       the case must not also meet another definition (focal and not generalized)
  bars        the code must be the person's most frequent diagnosis of that kind
  target      the control rule differs from the set's first definition
A faint dot means the rule does not apply. Two definitions with identical code rows can still differ
in the right block, which is exactly where some of them differ. Tiles show ICD-10 only; ICD-9 and
ICD-8 rules are in the endpoint file. Source colors pass the dataviz palette validator on white.

Sources: results/definitions/code_membership.tsv and endpoint_anatomy.tsv (src/02).
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

import viz_style as vs
from itoju_labels import ENDPOINT_LABEL

ROOT = vs.ROOT
HOSP = "#4b57c9"
PRIM = "#12998c"
SETS = [
    ("Epilepsy", ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"]),
    ("Sleep apnoea", ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"]),
    ("Insomnia", ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"]),
    ("Constipation", ["K11_CONSTIPATION", "K11_OTHFUNC"]),
    ("ADHD", ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"]),
]
RULES = [("drug", "drug purchase"), ("reimb", "reimbursement"), ("excl", "excludes another\ndefinition"),
         ("mode", "most frequent\ncode"), ("ctrl", "control rule\ndiffers")]
RULE_GAP = 1.5
TOP = 6.4          # header band for the rule names, drawn once
BLOCK_HEAD = 2.4   # per condition: its name and its rotated code labels
UNIT_IN = 0.14     # inches per tile


def glyph(ax, kind, x, y, on=True):
    """One rule glyph centered at (x, y), in tile units."""
    ink = vs.INK_2
    if not on:
        ax.add_patch(Circle((x, y), 0.07, facecolor=vs.GRID, edgecolor="none", zorder=3))
        return
    if kind == "drug":
        ax.add_patch(FancyBboxPatch((x - 0.36, y - 0.15), 0.72, 0.30, boxstyle="round,pad=0,rounding_size=0.15",
                                    facecolor="white", edgecolor=ink, lw=0.8, zorder=3))
        ax.add_patch(FancyBboxPatch((x - 0.36, y - 0.15), 0.36, 0.30, boxstyle="round,pad=0,rounding_size=0.15",
                                    facecolor=ink, edgecolor=ink, lw=0.8, zorder=4))
    elif kind == "reimb":
        ax.add_patch(Rectangle((x - 0.34, y - 0.24), 0.68, 0.48, facecolor="white", edgecolor=ink, lw=0.8, zorder=3))
        ax.add_patch(Rectangle((x - 0.34, y - 0.16), 0.68, 0.1, facecolor=ink, edgecolor="none", zorder=4))
    elif kind == "excl":
        ax.add_patch(Circle((x, y), 0.28, facecolor="white", edgecolor=ink, lw=0.9, zorder=3))
        ax.plot([x - 0.2, x + 0.2], [y - 0.2, y + 0.2], color=ink, lw=0.9, zorder=4)
    elif kind == "mode":
        for dx, h, fill in ((-0.22, 0.26, False), (0.0, 0.56, True), (0.22, 0.34, False)):
            ax.add_patch(Rectangle((x + dx - 0.08, y + 0.28 - h), 0.16, h, facecolor=ink if fill else "white",
                                   edgecolor=ink, lw=0.7, zorder=3))
    elif kind == "ctrl":
        ax.add_patch(Circle((x, y), 0.29, facecolor="white", edgecolor=ink, lw=0.8, zorder=3))
        ax.add_patch(Circle((x, y), 0.11, facecolor=ink, edgecolor="none", zorder=4))


def tile(ax, x, y, hosp, prim):
    if hosp:
        ax.add_patch(Rectangle((x - 0.42, y - 0.42), 0.84, 0.84, facecolor=HOSP, edgecolor="none", zorder=2))
        if prim:
            # y grows downward in these axes, so y - 0.42 is the tile's top edge: a top-right corner
            ax.add_patch(Polygon([(x + 0.42, y - 0.42), (x + 0.42, y + 0.06), (x - 0.06, y - 0.42)], closed=True,
                                 facecolor=PRIM, edgecolor="white", lw=0.5, zorder=3))
    else:
        ax.add_patch(Rectangle((x - 0.42, y - 0.42), 0.84, 0.84, facecolor="white", edgecolor=vs.GRID, lw=0.5,
                               zorder=2))


def main():
    vs.apply()
    mem = pd.read_csv(ROOT / "results" / "definitions" / "code_membership.tsv", sep="\t")
    ana = pd.read_csv(ROOT / "results" / "definitions" / "endpoint_anatomy.tsv", sep="\t").set_index("endpoint")

    layout, y = [], TOP
    for name, defs in SETS:
        codes = sorted(set(mem[mem.endpoint.isin(defs)].icd10))
        codes = [c for c in codes if not (len(c) == 3 and any(o.startswith(c) and len(o) > 3 for o in codes))] or codes
        layout.append((name, defs, codes, y))
        y += BLOCK_HEAD + len(defs) + 0.6
    total = y - 0.6 + 0.3
    rx0 = max(len(c) for _, _, c, _ in layout) + 1.2
    x_max = rx0 + (len(RULES) - 1) * RULE_GAP + 0.7
    x_min = -0.7

    left_in, key_in = 1.30, 1.55
    body_w = (x_max - x_min) * UNIT_IN
    fig_w = vs.DOUBLE
    fig_h = total * UNIT_IN + 0.08
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([left_in / fig_w, 0.04 / fig_h, body_w / fig_w, total * UNIT_IN / fig_h])
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(total, 0)
    ax.axis("off")

    for j, (_, lab) in enumerate(RULES):
        ax.text(rx0 + j * RULE_GAP, TOP - 0.5, lab, rotation=90, ha="center", va="bottom", fontsize=7,
                color=vs.INK_2, linespacing=1.05)
    ax.text(rx0 - 0.75, TOP - 0.5, "rules", ha="right", va="bottom", fontsize=7, color=vs.MUTED)
    ax.text(0.0 - 0.42, TOP - 0.5, "ICD-10 codes", ha="left", va="bottom", fontsize=7, color=vs.MUTED)

    for name, defs, codes, y0 in layout:
        ax.text(x_min, y0 + 0.9, name, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK)
        for k, c in enumerate(codes):
            ax.text(k, y0 + BLOCK_HEAD - 0.5, f"{c[:3]}.{c[3:]}" if len(c) > 3 else c, rotation=90, ha="center",
                    va="bottom", fontsize=7, color=vs.MUTED)
        ref_ctrl = str(ana.loc[defs[0], "control_exclude"]) + str(ana.loc[defs[0], "control_conditions"])
        for i, e in enumerate(defs):
            yy = y0 + BLOCK_HEAD + i
            a = ana.loc[e]
            ax.text(x_min, yy, ENDPOINT_LABEL[e], ha="right", va="center", fontsize=7, color=vs.INK_2)
            m = mem[mem.endpoint == e]
            hd = set(m[m.source.isin(["HD", "COD"])].icd10)
            op = set(m[m.source == "OUTPAT"].icd10)
            for k, c in enumerate(codes):
                tile(ax, k, yy, c in hd or any(h.startswith(c) for h in hd), c in op or any(o.startswith(c) for o in op))
            sources = str(a.sources).split()
            flags = {
                "drug": "KELA_ATC" in sources,
                "reimb": "KELA_REIMB" in sources,
                "excl": isinstance(a.case_conditions, str) and a.case_conditions.strip() != "",
                "mode": bool(a.mode_rule),
                "ctrl": (str(a.control_exclude) + str(a.control_conditions)) != ref_ctrl,
            }
            for j, (kind, _) in enumerate(RULES):
                glyph(ax, kind, rx0 + j * RULE_GAP, yy, flags[kind])

    # key, a vertical list on the right, every mark drawn by the same functions as the figure
    kx0 = (left_in + body_w + 0.25) / fig_w
    kw = 1 - kx0 - 0.01
    items = [("tile", (True, False), "code counted in hospital\nor death records"),
             ("tile", (True, True), "also counted from\nprimary care"),
             ("tile", (False, False), "code not counted"),
             *[("glyph", k, lab.replace("\n", " ")) for k, lab in RULES],
             ("off", None, "rule does not apply")]
    key_units = 1.6 + 1.55 * len(items)
    kax = fig.add_axes([kx0, 1 - (0.2 + key_units * UNIT_IN) / fig_h, kw, key_units * UNIT_IN / fig_h])
    kax.set_xlim(-0.6, kw * fig_w / UNIT_IN - 0.6)
    kax.set_ylim(key_units, 0)
    kax.axis("off")
    kax.text(-0.5, 0.5, "Key", fontsize=8, color=vs.INK, va="center")
    for i, (kind, arg, lab) in enumerate(items):
        yy = 1.9 + 1.55 * i
        if kind == "tile":
            tile(kax, 0, yy, *arg)
        elif kind == "glyph":
            glyph(kax, arg, 0, yy, True)
        else:
            glyph(kax, "none", 0, yy, False)
        kax.text(0.9, yy, lab, va="center", fontsize=7, color=vs.INK_2, linespacing=1.1)
    vs.save(fig, "fig02_anatomy")


if __name__ == "__main__":
    main()
