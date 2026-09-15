"""Figure 18. Sample overlap shows up where it should, in the LDSC intercept, and nowhere else.

The cross-trait LDSC intercept absorbs shared samples between two GWAS, which keeps rg unbiased.
Each row is one psychiatric GWAS; each dot is its intercept with one of the FinnGen endpoints, laid out
as a beeswarm so no dot hides another. The pale band around zero spans plus or minus 1.96 times the
median intercept SE: where the intercepts of a GWAS with no shared people should fall. Autism,
schizophrenia, bipolar disorder and clinically defined depression contain no FinnGen samples. The full
depression GWAS and its EHR subset include FinnGen R5, and the PTSD Freeze 3 cohort list includes
FinnGen ("fing"), so their intercepts should rise off zero. Black tick: median.

Source: results/ldsc/rg.tsv (standard two-step LDSC).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs
from itoju_labels import TRAIT_LABEL

ROOT = vs.ROOT
BLOCKS = [("No FinnGen samples", ["ASD", "SCZ", "BIP", "MDD_Clin"]),
          ("Includes FinnGen samples", ["PTSD", "MDD", "MDD_EHR"])]
SLEEVE = "#e6e1d6"
DOT_PT = 4.0
AX_LEFT, AX_W = 0.40, 0.57


def swarm(xs_pt, d):
    """Offsets (points) that keep dots of diameter d from touching, closest to the row line first."""
    placed, out = [], np.zeros(len(xs_pt))
    for i in np.argsort(xs_pt):
        x = xs_pt[i]
        for k in range(0, 80):
            off = ((k + 1) // 2) * d * (1 if k % 2 else -1) if k else 0.0
            if all((x - px) ** 2 + (off - py) ** 2 >= (d * 0.98) ** 2 for px, py in placed):
                break
        placed.append((x, off))
        out[i] = off
    return out


def main():
    vs.apply()
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    endpoints = set(rg[rg.trait1 == "ASD"].trait2) - {"SCZ", "BIP", "MDD", "PTSD", "MDD_EHR", "MDD_Clin"}
    traits = [t for _, ts in BLOCKS for t in ts]
    sub = rg[rg.trait1.isin(traits) & rg.trait2.isin(endpoints)].dropna(subset=["gcov_int"])
    band = 1.96 * float(sub.gcov_int_se.median())
    xlo = np.floor(sub.gcov_int.min() / 0.01) * 0.01 - 0.005
    xhi = np.ceil(sub.gcov_int.max() / 0.01) * 0.01 + 0.005

    # everything vertical is laid out in points, so each row gets exactly the height its swarm needs
    pt_per_x = AX_W * vs.SINGLE * 72 / (xhi - xlo)
    offs = {t: swarm(sub[sub.trait1 == t].gcov_int.values * pt_per_x, DOT_PT) for t in traits}
    y, heads, order = 8.0, [], []
    for block, ts in BLOCKS:
        heads.append((block, y))
        y += 12.0
        for t in ts:
            ext = np.abs(offs[t]).max() + DOT_PT / 2
            y += max(ext, 7.0)
            order.append((t, y, ext))
            y += max(ext, 7.0) + 5.0
        y += 6.0
    axes_pt, top_pt, bottom_pt = y, 34.0, 36.0
    fig_h = (axes_pt + top_pt + bottom_pt) / 72
    fig = plt.figure(figsize=(vs.SINGLE, fig_h))
    ax = fig.add_axes([AX_LEFT, bottom_pt / 72 / fig_h, AX_W, axes_pt / 72 / fig_h])
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(axes_pt, 0)
    ax.axvspan(-band, band, color=SLEEVE, lw=0, zorder=0)
    ax.axvline(0, color=vs.INK_2, lw=0.5, zorder=1)
    for block, yy in heads:
        ax.text(-0.02, yy, block, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.6,
                fontweight="bold", color=vs.INK)
    for t, yy, ext in order:
        s_ = sub[sub.trait1 == t]
        face = "white" if t in vs.HOLLOW else vs.TRAIT[t]
        ax.scatter(s_.gcov_int, yy + offs[t], s=DOT_PT ** 2 * 0.8, facecolor=face, edgecolor=vs.TRAIT[t],
                   linewidths=0.8, zorder=3)
        med = float(np.median(s_.gcov_int))
        ax.plot([med, med], [yy - min(ext, 9), yy + min(ext, 9)], color=vs.INK, lw=1.2, zorder=4)
        ax.text(-0.02, yy, TRAIT_LABEL[t], transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=7.5, color=vs.INK_2)
        print(f"  {t:9s} n {len(s_)} median {med:.4f} range {s_.gcov_int.min():.4f} to {s_.gcov_int.max():.4f}")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("cross-trait LDSC intercept", fontsize=7.5)
    fig.text(0.02, 1 - 6 / 72 / fig_h, "One dot per FinnGen endpoint.  Black tick: median.\n"
             f"Pale band: where no shared people should put it (1.96 x median SE, {band:.3f}).",
             fontsize=7, color=vs.INK_2, va="top", linespacing=1.35)
    vs.save(fig, "fig18_intercepts")


if __name__ == "__main__":
    main()
