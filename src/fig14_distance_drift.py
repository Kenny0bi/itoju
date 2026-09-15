"""Figure 14. Do definitions that differ more give genetic correlations that differ more?

Each point is one comparison between two definitions of the same condition, for one psychiatric
trait (135 comparisons). x is how different the two definitions are; y is how far rg moved.
  a  x = the pre-specified definition distance: 1 - mean(ICD-10 code Jaccard, data-source Jaccard,
       1 - any rule difference), fixed in the lab notebook before any rg was estimated.
  b  x = 1 - SNOMED CT hierarchy-aware Jaccard (concepts plus ancestors within two levels).
Spearman correlations are printed per trait. The comparisons share traits and definitions, so the
points are not independent and the correlations are descriptive, not tests.

Sources: results/definition_effect/pairwise_delta.tsv, results/omop/pairwise_omop.tsv.
Writes results/definition_effect/distance_vs_delta.tsv.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import viz_style as vs
from itoju_labels import TRAIT_LABEL

ROOT = vs.ROOT
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
FAMILY_SHAPE = {"Sleep apnoea": "o", "Insomnia": "D", "Constipation": "s", "ADHD": "^", "Epilepsy": "v"}


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

    stats = []
    fig, axes = plt.subplots(1, 2, figsize=(vs.DOUBLE, 3.1), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.78, bottom=0.15, wspace=0.08)
    for ax, (col, xlabel, letter) in zip(axes, [("distance", "definition distance (codes, sources, rules)", "a  Pre-specified distance"),
                                                 ("snomed_distance", "SNOMED CT distance (1 - hierarchy-aware Jaccard)", "b  Clinical-meaning distance (exploratory)")]):
        rng = np.random.default_rng(7)
        for r in pw.itertuples():
            jitter = rng.uniform(-0.018, 0.018)
            ax.scatter(getattr(r, col) + jitter, r.abs_delta, marker=FAMILY_SHAPE[r.family], s=20,
                       facecolor=vs.TRAIT[r.shared], edgecolor="white", linewidths=0.4,
                       alpha=0.6 if r.family == "Epilepsy" else 0.95, zorder=3)
        lines = []
        for t in TRAITS:
            s = pw[pw.shared == t]
            rho, p = spearmanr(s[col], s.abs_delta)
            stats.append({"distance": col, "shared": t, "n": len(s), "spearman_rho": rho, "p_descriptive": p})
            lines.append(f"{TRAIT_LABEL[t]}: {rho:+.2f}")
        rho_all, p_all = spearmanr(pw[col], pw.abs_delta)
        stats.append({"distance": col, "shared": "all", "n": len(pw), "spearman_rho": rho_all, "p_descriptive": p_all})
        tx, ha = (0.98, "right") if col == "distance" else (0.30, "left")
        ax.text(tx, 0.97, "Spearman rho\n" + "\n".join(lines) + f"\nall: {rho_all:+.2f}", transform=ax.transAxes,
                va="top", ha=ha, fontsize=7, color=vs.INK_2)
        ax.set_xlabel(xlabel)
        ax.set_xlim((0.28, 0.88) if col == "distance" else (-0.04, 1.02))
        ax.set_title(f"{letter}", loc="left", fontsize=8.5)
    axes[0].set_ylabel("|difference in rg|")
    axes[0].set_ylim(0, 0.5)
    pd.DataFrame(stats).to_csv(ROOT / "results" / "definition_effect" / "distance_vs_delta.tsv", sep="\t", index=False)

    def key_row(items, y):
        x = 0.08
        for lab, style in items:
            fig.add_artist(plt.Line2D([x], [y], ls="", transform=fig.transFigure, **style))
            fig.text(x + 0.009, y, lab, fontsize=7, color=vs.INK_2, va="center")
            x += 0.022 + 0.0105 * len(lab)
    key_row([(TRAIT_LABEL[t], dict(marker="o", markersize=4.5, markerfacecolor=vs.TRAIT[t], markeredgecolor=vs.TRAIT[t]))
             for t in TRAITS], 0.955)
    key_row([(f if f != "Epilepsy" else "Epilepsy (drawn lighter)",
              dict(marker=m, markersize=4.5, markerfacecolor="white", markeredgecolor=vs.INK_2))
             for f, m in FAMILY_SHAPE.items()], 0.895)
    vs.save(fig, "fig14_distance_drift")
    print(pd.DataFrame(stats).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
