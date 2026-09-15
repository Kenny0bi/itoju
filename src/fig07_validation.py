"""Figure 7. Validation: the pipeline lands inside published estimates before it is trusted.

Single column, two panels stacked.
a  Each published estimate is drawn the way every comparison in this study is drawn. The pale sleeve,
   centered on the published value, is the range that value allows: its 95% interval, or, when the
   paper prints no SE, the tolerance of plus or minus 0.05 fixed in the analysis plan. The stem runs
   from the published value to this pipeline's estimate, so a stem that stays inside its sleeve
   reproduces the paper.
   1. Autism SNP heritability, liability scale, K = 0.012 (Grove et al. 2019: 0.118, SE 0.010)
   2. Autism x ADHD genetic correlation (Grove et al. 2019: 0.360; no SE printed)
b  Does Finnish registry autism measure the same condition as the Danish autism GWAS? Each rg is
   drawn as a quantile dot plot of its sampling distribution: 20 dots, each one twentieth of the
   probability, against rg = 1, the value if both measured the same genetics.

Sources: results/ldsc/h2_ASD_Grove2019.log, rg_validation_ASD_ADHD.log, rg.tsv.
"""
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

import viz_style as vs

ROOT = vs.ROOT
Z = 1.96
SLEEVE = "#e6e1d6"
ASD = vs.TRAIT["ASD"]
BW = 0.06
NDOT = 20
FIG_H = 4.3


def grab(path, label):
    t = (ROOT / "results" / "ldsc" / path).read_text()
    m = re.search(re.escape(label) + r":\s*(-?[\d.]+)\s*\(([\d.]+)\)", t)
    return float(m.group(1)), float(m.group(2))


def main():
    vs.apply()
    h2, h2se = grab("h2_ASD_Grove2019.log", "Total Liability scale h2")
    rg, rgse = grab("rg_validation_ASD_ADHD.log", "Genetic Correlation")
    rows = [("Autism SNP heritability", "liability scale, K = 0.012", h2, 0.118, Z * 0.010, "95% interval"),
            ("Autism x ADHD", "genetic correlation", rg, 0.360, 0.05, "tolerance")]
    rgt = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    fin = [("FinnGen autism", "F84.0, F84.5; 888 cases", "KRA_PSY_AUTISM_EXMORE"),
           ("FinnGen autism spectrum", "F84; 1,179 cases", "KRA_PSY_DEVWIDE_EXMORE")]

    fig = plt.figure(figsize=(vs.SINGLE, FIG_H))
    fig.text(0.02, 0.985, "a  Reproducing the published autism results", fontsize=8.5, color=vs.INK, va="top")
    fig.text(0.02, 0.935, "Pale sleeve: the range the paper allows.  Dot: this pipeline.", fontsize=7,
             color=vs.INK_2, va="top")
    a = fig.add_axes([0.43, 0.64, 0.30, 0.22])
    for i, (lab, sub, est, pub, half, kind) in enumerate(rows):
        y = len(rows) - 1 - i
        d = est - pub
        a.hlines(y, -half, half, colors=SLEEVE, linewidth=8, zorder=1)
        a.hlines(y, 0, d, colors=ASD, linewidth=1.5, zorder=3)
        a.scatter(d, y, s=22, color=ASD, linewidths=0, zorder=4)
        a.text(-0.11, y + 0.14, lab, ha="right", va="center", fontsize=7.5, color=vs.INK)
        a.text(-0.11, y - 0.2, sub, ha="right", va="center", fontsize=7, color=vs.INK_2)
        a.text(0.11, y + 0.14, f"{est:.3f} vs {pub:.3f}", ha="left", va="center", fontsize=7, color=vs.INK)
        a.text(0.11, y - 0.2, kind, ha="left", va="center", fontsize=7, color=vs.INK_2)
    a.axvline(0, color=vs.INK_2, lw=0.5, zorder=2)
    a.set_xlim(-0.1, 0.1)
    a.set_ylim(-0.6, len(rows) - 0.4)
    a.set_yticks([])
    a.spines["left"].set_visible(False)
    a.set_xticks([-0.1, 0, 0.1])
    a.set_xticklabels(["-0.1", "0", "+0.1"])
    a.set_xlabel("this pipeline minus published", fontsize=7.5)

    fig.text(0.02, 0.475, "b  Finnish registry autism agrees with the autism GWAS", fontsize=8.5, color=vs.INK,
             va="top")
    fig.text(0.02, 0.425, "Each dot: 1/20 of the estimate's probability.", fontsize=7, color=vs.INK_2, va="top")
    x_lo, x_hi = 0.2, 1.5
    left_b, width_b, row_h = 0.43, 0.54, 0.115
    for i, (lab, sub, e) in enumerate(fin):
        r = rgt[(rgt.trait1 == "ASD") & (rgt.trait2 == e)].iloc[0]
        bottom = 0.265 - i * (row_h + 0.05)
        b = fig.add_axes([left_b, bottom, width_b, row_h])
        w_pt = width_b * vs.SINGLE * 72
        h_pt = row_h * FIG_H * 72
        unit = BW / (x_hi - x_lo) * w_pt
        q = norm.ppf((np.arange(NDOT) + 0.5) / NDOT, r.rg, r.se)
        bins = np.round(q / BW).astype(int)
        xs, hs = [], []
        for k in np.unique(bins):
            n = int((bins == k).sum())
            xs += [k * BW] * n
            hs += list(np.arange(n) + 0.5)
        b.scatter(xs, hs, s=(unit * 0.86) ** 2, color=ASD, linewidths=0, zorder=3)
        b.axvline(1.0, color=vs.INK_2, lw=0.7, ls=(0, (2, 2)), zorder=1)
        b.set_xlim(x_lo, x_hi)
        b.set_ylim(0, h_pt / unit)
        b.set_yticks([])
        b.spines["left"].set_visible(False)
        b.text(-0.04, 0.72, lab, transform=b.transAxes, ha="right", va="center", fontsize=7.5, color=vs.INK)
        b.text(-0.04, 0.22, f"{sub}\nrg {r.rg:.2f} (SE {r.se:.2f})", transform=b.transAxes, ha="right",
               va="center", fontsize=7, color=vs.INK_2, linespacing=1.15)
        if i == 0:
            b.text(1.02, h_pt / unit * 0.98, "same genetics", ha="left", va="top", fontsize=7, color=vs.INK_2)
            b.set_xticks([])
            b.spines["bottom"].set_color(vs.GRID)
        else:
            b.set_xticks([0.5, 1.0, 1.5])
            b.set_xticklabels(["0.5", "1", "1.5"])
            b.set_xlabel("rg with iPSYCH-PGC autism", fontsize=7.5)
    vs.save(fig, "fig07_validation")


if __name__ == "__main__":
    main()
