"""Figure 9. The genetic correlation quilt: seven psychiatric traits against 28 FinnGen endpoints.

Each cell is a dot. Its color is rg (diverging blue to red, capped at +/-1). Its area grows with
precision (1/SE), so an estimate with a wide standard error shrinks instead of shouting at the
same volume as a precise one. A thin ink ring marks p < 0.05. Rows are grouped by condition, with
the definitions of one condition stacked together, so a reader can scan down a column and see
whether switching definition moves the color.

Source: results/ldsc/rg.tsv (standard two-step LDSC, full SNP sets).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs
from itoju_labels import GROUPS

# depression and its two definitions sit together under one bracket
COLS = ["ASD", "SCZ", "BIP", "PTSD", "MDD", "MDD_EHR", "MDD_Clin"]
COL_LABEL = {"ASD": "ASD", "SCZ": "SCZ", "BIP": "BIP", "PTSD": "PTSD", "MDD": "all", "MDD_EHR": "EHR", "MDD_Clin": "clin."}

ROOT = vs.ROOT


def main():
    vs.apply()
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    rg = rg[rg.trait1.isin(COLS)]
    cell = {(r.trait1, r.trait2): r for r in rg.itertuples()}

    rows, headers = [], []
    y = 0.0
    for gname, items in GROUPS:
        headers.append((gname, y))
        y += 0.9
        for e, lab in items:
            rows.append((e, lab, y))
            y += 1.0
        y += 0.35
    height_units = y

    fig = plt.figure(figsize=(vs.SINGLE, 0.155 * height_units + 0.95))
    ax = fig.add_axes([0.42, 0.075, 0.56, 0.86])
    cmap = vs.cmap_div()
    CAP = 0.8  # every quilt value lies within -0.25 to 0.80, so +/-1 would wash the colors out
    norm = plt.Normalize(-CAP, CAP)
    xs = {t: i for i, t in enumerate(COLS)}

    max_area = 62.0
    for e, lab, yy in rows:
        for t in COLS:
            r = cell.get((t, e))
            if r is None or pd.isna(r.rg):
                ax.text(xs[t], yy, "NA", ha="center", va="center", fontsize=7, color=vs.MUTED)
                continue
            prec = min(1.0, 0.05 / r.se)  # SE 0.05 or smaller draws at full size
            area = max_area * prec
            ax.scatter(xs[t], yy, s=area, c=[cmap(norm(np.clip(r.rg, -CAP, CAP)))],
                       edgecolors=vs.INK_2 if r.p < 0.05 else "none", linewidths=0.4, zorder=3)
        ax.text(-0.75, yy, lab, ha="right", va="center", fontsize=7.2, color=vs.INK_2)

    for gname, yy in headers:
        ax.text(-0.75, yy, gname, ha="right", va="center", fontsize=7.6, color=vs.INK, fontweight="bold")
        ax.plot([-0.5, len(COLS) - 0.5], [yy + 0.45, yy + 0.45], color=vs.GRID, lw=0.5, zorder=1)

    # column headers: a small trait-colored dot beside ink text (text never wears the trait color)
    for t, x in xs.items():
        ax.text(x, -1.25, COL_LABEL[t], ha="center", va="bottom", fontsize=7, color=vs.INK_2)
        face = "white" if t in vs.HOLLOW else vs.TRAIT[t]
        ax.scatter(x, -0.55, s=14, facecolor=face, edgecolor=vs.TRAIT[t], linewidths=1.0, clip_on=False, zorder=4)
    x0, x1 = xs["MDD"] - 0.4, xs["MDD_Clin"] + 0.4
    ax.plot([x0, x0, x1, x1], [-2.35, -2.55, -2.55, -2.35], color=vs.MUTED, lw=0.6, clip_on=False)
    ax.text((x0 + x1) / 2, -2.75, "Depression", ha="center", va="bottom", fontsize=7, color=vs.INK_2)

    ax.set_xlim(-0.5, len(COLS) - 0.5)
    ax.set_ylim(height_units, -3.4)
    ax.axis("off")

    # color key
    cax = fig.add_axes([0.10, 0.028, 0.34, 0.010])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-0.8, -0.4, 0, 0.4, 0.8])
    cb.set_ticklabels(["-0.8", "-0.4", "0", "0.4", "0.8+"])
    cb.ax.tick_params(labelsize=7, length=2, colors=vs.MUTED)
    cb.outline.set_visible(False)
    cax.set_title("genetic correlation, rg", fontsize=7, color=vs.INK_2, pad=2)

    # size and ring key, drawn as real marks in a small key panel (not a legend of stand-in swatches)
    kax = fig.add_axes([0.50, 0.008, 0.48, 0.04])
    kax.set_xlim(0, 5.0)
    kax.set_ylim(0, 1)
    kax.axis("off")
    for i, (f, lab) in enumerate(((1.0, "SE 0.05"), (0.5, "0.10"), (0.2, "0.25"))):
        x = [0.15, 1.55, 2.45][i]
        kax.scatter(x, 0.5, s=max_area * f, c=[cmap(norm(0.5))], edgecolors="none")
        kax.text(x + 0.2, 0.5, lab, va="center", fontsize=7, color=vs.INK_2)
    kax.scatter(3.4, 0.5, s=max_area * 0.5, c=[cmap(norm(0.5))], edgecolors=vs.INK_2, linewidths=0.4)
    kax.text(3.6, 0.5, "p < 0.05", va="center", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig09_rg_quilt")


if __name__ == "__main__":
    main()
