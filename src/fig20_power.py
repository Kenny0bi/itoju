"""Figure 20. What buys the power to see a definition effect: shared cases.

The smallest detectable difference between two genetic correlations is
    MDD = (z_0.975 + z_0.80) * sqrt(se1^2 + se2^2 - 2 * rho * se1 * se2),
where rho is the correlation between the two estimates' sampling errors. Two definitions that
share most of their cases have errors that move together (rho near 1), so most of the noise cancels.
Each point is one comparison between definitions of the same condition (135): x is its error
correlation from the block jackknife, y its minimum detectable difference (log scale). Gray curves
show the formula for two estimates with equal SE of 0.05, 0.10 and 0.20.

Source: results/definition_effect/pairwise_delta.tsv and rg_common.tsv.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

import viz_style as vs
from itoju_labels import TRAIT_LABEL

ROOT = vs.ROOT
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
FAMILY_SHAPE = {"Sleep apnoea": "o", "Insomnia": "D", "Constipation": "s", "ADHD": "^", "Epilepsy": "v"}
K = norm.ppf(0.975) + norm.ppf(0.80)


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw = pw[pw.family != "MDD definition (positive control)"]
    fig, ax = plt.subplots(figsize=(vs.SINGLE, 3.1))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.875, bottom=0.14)
    rho = np.linspace(-0.2, 0.995, 200)
    for se in (0.05, 0.10, 0.20):
        ax.plot(rho, K * se * np.sqrt(2 - 2 * rho), color=vs.GRID, lw=1.2, zorder=1)
    for r in pw.itertuples():
        ax.scatter(r.error_corr, r.min_detectable_delta, marker=FAMILY_SHAPE[r.family], s=20,
                   facecolor=vs.TRAIT[r.shared], edgecolor="white", linewidths=0.4,
                   alpha=0.6 if r.family == "Epilepsy" else 0.95, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 1.5)
    ax.set_yticks([0.01, 0.03, 0.1, 0.3, 1.0])
    ax.set_yticklabels(["0.01", "0.03", "0.1", "0.3", "1"])
    ax.set_xlim(-0.2, 1.0)
    ax.set_xlabel("correlation between the two estimates' errors")
    ax.set_ylabel("smallest detectable difference in rg")

    def key_row(items, y, x0=0.03):
        # spacing from label length at 7 pt on a 3.5 in wide figure (about 0.0155 figure widths per character)
        x = x0
        for lab, style in items:
            fig.add_artist(plt.Line2D([x], [y], ls="", transform=fig.transFigure, **style))
            fig.text(x + 0.018, y, lab, fontsize=7, color=vs.INK_2, va="center")
            x += 0.065 + 0.0165 * len(lab)
    key_row([(TRAIT_LABEL[t] if t != "SCZ" else "Schizo.", dict(marker="o", markersize=4, markerfacecolor=vs.TRAIT[t],
             markeredgecolor=vs.TRAIT[t])) for t in TRAITS], 0.97)
    key_row([({"Sleep apnoea": "Apnoea", "Constipation": "Constip.", "Epilepsy": "Epilepsy, lighter"}.get(f, f), dict(marker=m, markersize=4, markerfacecolor="white",
             markeredgecolor=vs.INK_2)) for f, m in FAMILY_SHAPE.items()], 0.925)
    ax.text(0.02, 0.03, "gray curves: two estimates with equal SE\n0.05, 0.10, 0.20 (bottom to top)", transform=ax.transAxes,
            fontsize=7, color=vs.INK_2, va="bottom")
    vs.save(fig, "fig20_power")


if __name__ == "__main__":
    main()
