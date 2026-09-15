"""Figure 9. The genetic correlation quilt, with its headline pattern read first.

Top: the takeaway in plain words, and a headline strip: for each condition, the median genetic
correlation across its definitions with each psychiatric trait (one dot per trait, in the trait's color).
The sentence is built from those medians and checked in code.

Below: the full matrix. Each tile is one genetic correlation between a psychiatric GWAS (column) and a
FinnGen endpoint (row). Color follows a value-suppressing uncertainty palette: the more precise an
estimate, the finer its color steps (8 steps when SE is at most 0.05); as the SE grows, neighbouring
values merge (4, then 2 steps) and the color fades toward gray, and above SE 0.20 every value is the
same gray. The fan key reads outward from uncertain (center) to precise (outer ring) and left to right
from negative to positive rg.

Source: results/ldsc/rg.tsv (standard two-step LDSC, full SNP sets).
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, Wedge

import viz_style as vs
from itoju_labels import GROUPS, TRAIT_SHORT

ROOT = vs.ROOT
COLS = ["ASD", "SCZ", "BIP", "PTSD", "MDD", "MDD_EHR", "MDD_Clin"]
FIVE = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
COL_LABEL = {"ASD": "ASD", "SCZ": "SCZ", "BIP": "BIP", "PTSD": "PTSD", "MDD": "all", "MDD_EHR": "EHR",
             "MDD_Clin": "clin."}
CAP = 0.8
SE_EDGES = [0.05, 0.10, 0.20]          # level 0: SE <= 0.05, 1: <= 0.10, 2: <= 0.20, 3: above
SUPPRESS = [0.0, 0.2, 0.45, 0.88]      # how far each level's colors are pulled toward gray
GRAYISH = np.array(mpl.colors.to_rgb("#d9d6cf"))
RAMP = mpl.colors.LinearSegmentedColormap.from_list("itoju_quilt", [
    (0.0, "#1f58b5"), (0.25, "#8db1e6"), (0.5, "#f3f0ea"), (0.625, "#f8c7ab"),
    (0.75, "#ee8a62"), (0.875, "#d44a37"), (1.0, "#971b2c")])
PITCH = 0.12                           # inches per tile row
COLP = 0.235                           # inches per tile column
LEFT_IN = 1.62
RINGS = [(1.10, 1.38), (0.84, 1.10), (0.58, 0.84), (0.32, 0.58)]
A0, A1 = 158.0, 22.0
HEAD_ROW = 0.13                        # inches per headline row


def level(se):
    return int(np.searchsorted(SE_EDGES, se, side="left"))


def vsup(v, lev, cmap):
    nb = 8 >> lev
    i = min(nb - 1, int((np.clip(v, -CAP, CAP) + CAP) / (2 * CAP) * nb))
    vc = -CAP + (i + 0.5) * 2 * CAP / nb
    base = np.array(cmap((vc + CAP) / (2 * CAP))[:3])
    return tuple(base * (1 - SUPPRESS[lev]) + GRAYISH * SUPPRESS[lev])


def main():
    vs.apply()
    cmap = RAMP
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    rg = rg[rg.trait1.isin(COLS)]
    cell = {(r.trait1, r.trait2): r for r in rg.itertuples()}

    # headline numbers and the sentence they support
    conds = [(g, [e for e, _ in items]) for g, items in GROUPS if g != "Neighbouring conditions"]
    med = {g: rg[rg.trait1.isin(FIVE) & rg.trait2.isin(es)].groupby("trait1").rg.median() for g, es in conds}
    top = {t: max(conds, key=lambda c: med[c[0]][t])[0] for t in FIVE}
    low = {t: min(conds, key=lambda c: med[c[0]][t])[0] for t in FIVE}
    n_top = sum(top[t] in ("ADHD", "Insomnia") for t in FIVE)
    n_low = sum(low[t] in ("Epilepsy", "Sleep apnoea") for t in FIVE)
    epi = dict(conds)["Epilepsy"]
    asd_epi = rg[(rg.trait1 == "ASD") & rg.trait2.isin(epi)]
    print("  strongest:", top, "\n  weakest:", low, f"\n  autism x epilepsy max rg {asd_epi.rg.max():.3f} over {len(asd_epi)}")
    assert n_top == 4 and n_low == 5 and asd_epi.rg.max() < 0 and len(asd_epi) == 7
    takeaway = (f"ADHD and insomnia share the most genetics with {n_top} of 5 psychiatric\n"
                f"traits, epilepsy or sleep apnoea the least with all {n_low}. Autism and\n"
                f"epilepsy correlate below zero under all {len(asd_epi)} epilepsy definitions.")

    rows, headers, y = [], [], 0.0
    for gname, items in GROUPS:
        headers.append((gname, y))
        y += 0.95
        for e, lab in items:
            rows.append((e, lab, y))
            y += 1.0
        y += 0.3
    units = y - 0.3
    head_units = 3.0
    strip_top, strip_h = 0.74, len(conds) * HEAD_ROW
    matrix_top = strip_top + strip_h + 0.36
    bottom_in = 0.95
    fig_h = matrix_top + (units + head_units) * PITCH + bottom_in
    fig = plt.figure(figsize=(vs.SINGLE, fig_h))
    frac = lambda top_in, h_in: (fig_h - top_in - h_in) / fig_h

    fig.text(0.02, 1 - 0.04 / fig_h, takeaway, fontsize=7.5, color=vs.INK, va="top", linespacing=1.25)
    # trait key for the strip
    kx = 0.36
    for t in FIVE:
        fig.add_artist(plt.Line2D([kx], [frac(0.58, 0)], marker="o", markersize=4.2, color=vs.TRAIT[t], ls="",
                                  transform=fig.transFigure))
        fig.text(kx + 0.017, frac(0.58, 0), TRAIT_SHORT[t], fontsize=7, color=vs.INK_2, va="center")
        kx += 0.05 + 0.022 * len(TRAIT_SHORT[t])
    fig.text(0.02, frac(0.58, 0), "median rg per trait", fontsize=7, color=vs.INK_2, va="center")

    s = fig.add_axes([LEFT_IN / vs.SINGLE, frac(strip_top, strip_h), (vs.SINGLE - LEFT_IN - 0.08) / vs.SINGLE,
                      strip_h / fig_h])
    s.set_xlim(-0.25, 0.85)
    s.set_ylim(len(conds) - 0.5, -0.5)
    for i, (g, _) in enumerate(conds):
        s.plot([-0.25, 0.85], [i, i], color=vs.GRID, lw=0.5, zorder=0)
        s.text(-0.04, i, g, transform=s.get_yaxis_transform(), ha="right", va="center", fontsize=7, color=vs.INK_2)
        for t in FIVE:
            s.scatter(med[g][t], i, s=17, color=vs.TRAIT[t], edgecolor="white", linewidths=0.4, zorder=3)
    s.axvline(0, color=vs.MUTED, lw=0.5, zorder=1)
    s.set_yticks([])
    s.spines["left"].set_visible(False)
    s.set_xticks([0, 0.4, 0.8])
    s.set_xticklabels(["0", "0.4", "0.8"])
    s.tick_params(axis="x", labelsize=7, pad=1)
    s.set_xlabel("median rg across the condition's definitions", fontsize=7, labelpad=1)

    ax = fig.add_axes([LEFT_IN / vs.SINGLE, bottom_in / fig_h, len(COLS) * COLP / vs.SINGLE,
                       (units + head_units) * PITCH / fig_h])
    ax.set_xlim(0, len(COLS))
    ax.set_ylim(units, -head_units)
    ax.axis("off")
    counts = [0, 0, 0, 0]
    for e, lab, yy in rows:
        ax.text(-0.3, yy + 0.5, lab, ha="right", va="center", fontsize=7, color=vs.INK_2)
        missing = 0
        for j, t in enumerate(COLS):
            r = cell.get((t, e))
            if r is None or pd.isna(r.rg):
                missing += 1
                continue
            lev = level(r.se)
            counts[lev] += 1
            ax.add_patch(Rectangle((j + 0.06, yy + 0.06), 0.88, 0.88, facecolor=vsup(r.rg, lev, cmap),
                                   edgecolor="none", zorder=2))
        if missing == len(COLS):
            ax.text(len(COLS) / 2, yy + 0.5, "no estimate: h2 below 0", ha="center", va="center",
                    fontsize=7, color=vs.MUTED)
    print("  tiles per SE level:", counts)
    for gname, yy in headers:
        ax.text(-0.3, yy + 0.48, gname, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK)
    for j, t in enumerate(COLS):
        ax.text(j + 0.5, -0.95, COL_LABEL[t], ha="center", va="bottom", fontsize=7, color=vs.INK_2)
        face = "white" if t in vs.HOLLOW else vs.TRAIT[t]
        ax.scatter(j + 0.5, -0.4, s=14, facecolor=face, edgecolor=vs.TRAIT[t], linewidths=1.0, zorder=4)
    x0, x1 = COLS.index("MDD") + 0.1, COLS.index("MDD_Clin") + 0.9
    ax.plot([x0, x0, x1, x1], [-1.95, -2.15, -2.15, -1.95], color=vs.MUTED, lw=0.6)
    ax.text((x0 + x1) / 2, -2.3, "Depression", ha="center", va="bottom", fontsize=7, color=vs.INK_2)
    fig.text(0.02, frac(matrix_top + 0.12, 0), "Full matrix", fontsize=7.5, color=vs.INK, va="center")

    kax = fig.add_axes([0.015, 0.006, 0.52, (bottom_in - 0.25) / fig_h])
    kax.set_aspect("equal")
    kax.axis("off")
    kax.set_xlim(-1.62, 1.62)
    kax.set_ylim(-0.12, 1.72)
    for lev, (r_in, r_out) in enumerate(RINGS):
        nb = 8 >> lev
        for i in range(nb):
            a_hi = A0 - i * (A0 - A1) / nb
            a_lo = A0 - (i + 1) * (A0 - A1) / nb
            vc = -CAP + (i + 0.5) * 2 * CAP / nb
            kax.add_patch(Wedge((0, 0), r_out, a_lo, a_hi, width=r_out - r_in, facecolor=vsup(vc, lev, cmap),
                                edgecolor="white", lw=0.8))
    for v, lab, ha, va in ((-0.8, "-0.8", "right", "center"), (0.0, "0", "center", "bottom"),
                           (0.8, "0.8+", "left", "center")):
        a = np.deg2rad(A0 - (v + CAP) / (2 * CAP) * (A0 - A1))
        kax.text(1.46 * np.cos(a), 1.46 * np.sin(a), lab, ha=ha, va=va, fontsize=7, color=vs.INK_2)
    fig.text(0.545, (bottom_in - 0.3) / fig_h, "Hue: rg, blue negative, red positive.\n"
             "Depth of color: precision.\n"
             "Outer ring, SE to 0.05: 8 steps\n"
             "SE to 0.10: 4 steps; to 0.20: 2\n"
             "Center, SE above 0.20: gray", fontsize=7, color=vs.INK_2, va="top", linespacing=1.3)
    vs.save(fig, "fig09_rg_quilt")


if __name__ == "__main__":
    main()
