"""Figure 18. Sample overlap shows up where it should, in the LDSC intercept, and nowhere else.

The cross-trait LDSC intercept absorbs shared samples between two GWAS, which keeps rg unbiased.
Each row is one psychiatric GWAS; each dot is its intercept with one of the 28 FinnGen endpoints,
with a 95% interval. Autism, schizophrenia, bipolar disorder and clinically defined depression
contain no FinnGen samples, so their intercepts should sit on zero.
The full depression GWAS and its EHR subset include FinnGen R5, and the PTSD Freeze 3 cohort list
includes FinnGen ("fing"), so their intercepts should sit above zero.

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
Z = 1.96


def main():
    vs.apply()
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    endpoints = set(rg[rg.trait1 == "ASD"].trait2) - {"SCZ", "BIP", "MDD", "PTSD", "MDD_EHR", "MDD_Clin"}
    fig, ax = plt.subplots(figsize=(vs.SINGLE, 3.0))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.90, bottom=0.15)
    rng = np.random.default_rng(3)
    y, order = 0.0, []
    for block, traits in BLOCKS:
        ax.text(-0.02, y, block, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.6,
                fontweight="bold", color=vs.INK)
        y += 1.0
        for t in traits:
            order.append((t, y))
            y += 1.0
        y += 0.4
    for t, yy in order:
        s_ = rg[(rg.trait1 == t) & rg.trait2.isin(endpoints)].dropna(subset=["gcov_int"])
        jit = rng.uniform(-0.25, 0.25, len(s_))
        for r, j in zip(s_.itertuples(), jit):
            ax.plot([r.gcov_int - Z * r.gcov_int_se, r.gcov_int + Z * r.gcov_int_se], [yy + j, yy + j],
                    color=vs.TRAIT[t], lw=0.5, alpha=0.35, zorder=2)
        face = "white" if t in vs.HOLLOW else vs.TRAIT[t]
        ax.scatter(s_.gcov_int, yy + jit, s=11, facecolor=face, edgecolor=vs.TRAIT[t], linewidths=0.7, zorder=3)
        med = float(np.median(s_.gcov_int))
        ax.plot([med, med], [yy - 0.38, yy + 0.38], color=vs.INK, lw=1.0, zorder=4)
        ax.text(-0.02, yy, TRAIT_LABEL[t], transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=7.5, color=vs.INK_2)
    ax.axvline(0, color=vs.MUTED, lw=0.6, zorder=1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_ylim(y - 0.2, -0.7)
    ax.set_xlabel("cross-trait LDSC intercept (95% interval)")
    fig.text(0.02, 0.955, "One dot per FinnGen endpoint; black tick: median.", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig18_intercepts")
    for t in [t for _, ts in BLOCKS for t in ts]:
        s = rg[(rg.trait1 == t) & rg.trait2.isin(endpoints)].dropna(subset=["gcov_int"])
        print(f"{t:9s} n {len(s)} median intercept {s.gcov_int.median():.4f} range {s.gcov_int.min():.4f} to {s.gcov_int.max():.4f}")


if __name__ == "__main__":
    main()
