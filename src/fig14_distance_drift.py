"""Figure 14. Do definitions that differ more give genetic correlations that differ more?

Single column, two panels stacked. Each point is one comparison between two definitions of the same
condition, for one psychiatric trait (135 comparisons), in the color of that trait; epilepsy's 105 are
hollow. x is how different the two definitions are; y is how far rg moved.
  a  x = the pre-specified definition distance: 1 - mean(ICD-10 code Jaccard, data-source Jaccard,
       1 - any rule difference), fixed in the lab notebook before any rg was estimated.
  b  x = 1 - SNOMED CT hierarchy-aware Jaccard (concepts plus ancestors within two levels).
Definition pairs share exact distances, so comparisons fall in clusters; for each cluster of five or more,
a dark bar marks the median |difference| and a pale box the middle half. Above each panel, a strip places
the Spearman correlation for each trait (colored dots) and for all comparisons (dark bar) on one scale.
The comparisons share traits and definitions, so the points are not independent and the correlations
are descriptive, not tests.

Sources: results/definition_effect/pairwise_delta.tsv, results/omop/pairwise_omop.tsv.
Writes results/definition_effect/distance_vs_delta.tsv.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr

import viz_style as vs
from itoju_labels import TRAIT_SHORT

ROOT = vs.ROOT
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
BAND = "#e6e1d6"
W, H = vs.SINGLE, 6.62
LEFT, RIGHT = 0.52, 0.08
PANELS = [
    ("distance", "definition distance (codes, sources, rules)", "a  Codes, sources and rules (pre-specified)",
     (0.28, 0.88)),
    ("snomed_distance", "SNOMED CT distance (1 - hierarchy-aware Jaccard)", "b  Clinical meaning (exploratory)",
     (-0.04, 1.02)),
]


def box(top, h):
    return [LEFT / W, (H - top - h) / H, (W - LEFT - RIGHT) / W, h / H]


def y_of(top):
    return 1 - top / H


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw = pw[pw.family != "MDD definition (positive control)"].copy()
    om = pd.read_csv(ROOT / "results" / "omop" / "pairwise_omop.tsv", sep="\t")
    om = pd.concat([om[["def1", "def2", "snomed_hier_jaccard"]],
                    om.rename(columns={"def1": "def2", "def2": "def1"})[["def1", "def2", "snomed_hier_jaccard"]]])
    pw = pw.merge(om.drop_duplicates(["def1", "def2"]), on=["def1", "def2"], how="left")
    pw["snomed_distance"] = 1 - pw.snomed_hier_jaccard
    pw["abs_delta"] = pw.delta.abs()

    fig = plt.figure(figsize=(W, H))
    x = 0.06
    for t in TRAITS:
        fig.add_artist(plt.Line2D([x / W], [y_of(0.12)], marker="o", markersize=4.5, color=vs.TRAIT[t], ls="",
                                  transform=fig.transFigure))
        fig.text((x + 0.07) / W, y_of(0.12), TRAIT_SHORT[t], fontsize=7, color=vs.INK_2, va="center")
        x += 0.2 + 0.06 * len(TRAIT_SHORT[t])
    fig.text(0.06 / W, y_of(0.30), "Filled: sleep, constipation and ADHD.  Hollow: epilepsy.", fontsize=7,
             color=vs.INK_2, va="center")
    fig.text(0.06 / W, y_of(0.46), "Dark bar: median of each cluster of 5 or more; box: middle half.", fontsize=7,
             color=vs.INK_2, va="center")

    stats = []
    epi = (pw.family == "Epilepsy").values
    cols = np.array([vs.TRAIT[t] for t in pw.shared])
    for k, (col, xlabel, head, xlim) in enumerate(PANELS):
        t0 = 0.72 + k * 3.0
        rho_all, p_all = spearmanr(pw[col], pw.abs_delta)
        fig.text(0.06 / W, y_of(t0), head, fontsize=8.5, color=vs.INK, va="top")
        fig.text(0.06 / W, y_of(t0 + 0.2), f"Spearman rho: dots by trait, bar for all 135 ({rho_all:.2f})", fontsize=7,
                 color=vs.INK_2, va="top")

        sx = fig.add_axes(box(t0 + 0.44, 0.16))
        sx.set_xlim(-0.3, 1.0)
        sx.set_ylim(0, 1)
        for s in ("left", "right", "top"):
            sx.spines[s].set_visible(False)
        sx.set_yticks([])
        sx.set_xticks([0, 0.5, 1.0])
        sx.set_xticklabels(["0", "0.5", "1"])
        sx.tick_params(axis="x", labelsize=7, length=2, pad=1)
        for t in TRAITS:
            s = pw[pw.shared == t]
            rho, p = spearmanr(s[col], s.abs_delta)
            stats.append({"distance": col, "shared": t, "n": len(s), "spearman_rho": rho, "p_descriptive": p})
            sx.scatter(rho, 0.45, s=22, color=vs.TRAIT[t], edgecolor="white", linewidths=0.5, zorder=3)
        stats.append({"distance": col, "shared": "all", "n": len(pw), "spearman_rho": rho_all, "p_descriptive": p_all})
        sx.plot([rho_all, rho_all], [0.0, 0.95], color=vs.INK, lw=2.0, zorder=4)

        ax = fig.add_axes(box(t0 + 0.9, 1.68))
        jit = np.random.default_rng(7).uniform(-0.012, 0.012, len(pw))
        xv, yv = pw[col].values + jit, pw.abs_delta.values
        ax.scatter(xv[~epi], yv[~epi], s=14, c=list(cols[~epi]), edgecolor="white", linewidths=0.3, zorder=3)
        ax.scatter(xv[epi], yv[epi], s=11, facecolor="white", edgecolor=list(cols[epi]), linewidths=0.7, zorder=2)
        clusters = np.round(pw[col].values / 0.04)
        for c_ in np.unique(clusters):
            b = pw[clusters == c_]
            if len(b) < 5:
                continue
            lo, hi = b[col].min() - 0.016, b[col].max() + 0.016
            q1, q3 = b.abs_delta.quantile(0.25), b.abs_delta.quantile(0.75)
            ax.add_patch(Rectangle((lo, q1), hi - lo, q3 - q1, facecolor=BAND, edgecolor="none", zorder=1))
            ax.plot([lo, hi], [b.abs_delta.median()] * 2, color=vs.INK, lw=2.2, solid_capstyle="butt", zorder=4)
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 0.5)
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_xlabel(xlabel, fontsize=7.5, labelpad=2)
        ax.set_ylabel("|difference in rg|", fontsize=7.5)
    pd.DataFrame(stats).to_csv(ROOT / "results" / "definition_effect" / "distance_vs_delta.tsv", sep="\t", index=False)
    vs.save(fig, "fig14_distance_drift")
    print(pd.DataFrame(stats).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
