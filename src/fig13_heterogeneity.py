"""Figure 13. Where does the definition matter? Each p-value drawn as the tail it is.

Rows are the definition sets, columns the five psychiatric traits. Each cell asks whether rg with
that trait varies across the definitions of that condition more than the jackknife covariance allows
(Cochran's Q with the full covariance matrix, df = number of definitions - 1). The curve is the
chi-square distribution Q would follow if the definitions made no difference. The black tick is the
observed Q, and the shaded tail beyond it is the p-value, drawn as the area it is. Tails with p < 0.05
wear the trait's color. A Q beyond the drawn range is shown by an arrow off the right edge.

Source: results/definition_effect/heterogeneity.tsv
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as mc
from scipy.stats import chi2

import viz_style as vs
from itoju_labels import TRAIT_SHORT

ROOT = vs.ROOT
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
SETS = [("Sleep apnoea", 3), ("Insomnia", 2), ("Constipation", 2), ("ADHD", 2), ("Epilepsy", 7)]
GRAY = "#d3cfc6"


def main():
    vs.apply()
    h = pd.read_csv(ROOT / "results" / "definition_effect" / "heterogeneity.tsv", sep="\t")
    cell = {(r.shared, r.family): r for r in h.itertuples()}

    fig = plt.figure(figsize=(vs.SINGLE, 3.55))
    left, right, top, bottom = 0.235, 0.995, 0.855, 0.075
    cw, ch = (right - left) / len(TRAITS), (top - bottom) / len(SETS)
    for i, (fam, k) in enumerate(SETS):
        df = k - 1
        xmax = chi2.ppf(0.995, df) * 1.15
        x = np.linspace(0, xmax, 400)[1:]
        dens = chi2.pdf(x, df)
        ymax = chi2.pdf(0.35, df) * 1.25 if df == 1 else dens.max() * 1.25
        yc = top - (i + 0.5) * ch
        fig.text(left - 0.015, yc + 0.012, fam, ha="right", va="center", fontsize=7.5, color=vs.INK)
        fig.text(left - 0.015, yc - 0.03, f"{k} definitions, df {df}", ha="right", va="center", fontsize=7,
                 color=vs.MUTED)
        for j, t in enumerate(TRAITS):
            r = cell[(t, fam)]
            ax = fig.add_axes([left + j * cw + 0.006, top - (i + 1) * ch + 0.012, cw - 0.012, ch - 0.066])
            ax.plot(x, np.minimum(dens, ymax), color=vs.INK_2, lw=0.6, zorder=2)
            clear = r.p < 0.05
            # the tail is the trait's color lightened toward white, so the tick stays the strongest mark
            col = tuple(0.45 * np.array(mc.to_rgb(vs.TRAIT[t])) + 0.55)
            tick = vs.TRAIT[t]
            if r.Q < xmax:
                m = x >= r.Q
                ax.fill_between(x[m], 0, np.minimum(dens[m], ymax), color=col, lw=0, zorder=1)
                ax.plot([r.Q, r.Q], [0, ymax * 0.8], color=tick, lw=1.4 if clear else 0.9, zorder=3)
                if clear:
                    # a significant tail is a sliver too thin to see, so the tick itself carries the color
                    ax.scatter(r.Q, ymax * 0.8, s=14, color=tick, linewidths=0, zorder=4, clip_on=False)
            else:
                ax.annotate("", xy=(xmax, ymax * 0.45), xytext=(xmax * 0.62, ymax * 0.45),
                            arrowprops=dict(arrowstyle="-|>", color=tick, lw=1.2, mutation_scale=6), zorder=3)
            ax.axhline(0, color=vs.MUTED, lw=0.5, zorder=0)
            ptxt = "p < 0.001" if r.p < 0.001 else f"p {r.p:.3f}"
            ax.text(0.5, -0.06, ptxt, transform=ax.transAxes, ha="center", va="top", fontsize=7,
                    color=vs.INK if clear else vs.INK_2, zorder=4)
            ax.set_xlim(0, xmax)
            ax.set_ylim(0, ymax * 1.02)
            ax.axis("off")
    for j, t in enumerate(TRAITS):
        xc = left + (j + 0.5) * cw
        fig.text(xc, top + 0.04, TRAIT_SHORT[t], ha="center", va="center", fontsize=7.5, color=vs.INK_2)
        fig.add_artist(plt.Line2D([xc], [top + 0.075], marker="o", markersize=4, color=vs.TRAIT[t], ls="",
                                  transform=fig.transFigure))
    fig.text(0.02, 0.975, "Curve: Q if the definitions made no difference.  Tick: observed Q.", fontsize=7,
             color=vs.INK_2)
    fig.text(0.02, 0.02, "Shaded tail: the p-value.  Dot on the tick: p < 0.05.  Arrow: Q off scale.",
             fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig13_heterogeneity")


if __name__ == "__main__":
    main()
