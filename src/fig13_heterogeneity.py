"""Figure 13. Where does the definition matter? Heterogeneity of rg across definitions.

Rows are the definition sets, columns the five psychiatric traits. Each dot summarises one test:
does rg with that trait vary across the definitions of that condition more than the jackknife
covariance allows (Cochran's Q with the full covariance matrix, df = number of definitions - 1)?
  fill   -log10 p of Q, one blue ramp, capped at 8 ("8+")
  size   the range of rg across the definitions (largest minus smallest)
  ring   p < 0.05
Printed beside the dot: Q and df for tests with p < 0.05.

Source: results/definition_effect/heterogeneity.tsv
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs
from itoju_labels import TRAIT_SHORT

ROOT = vs.ROOT
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
SETS = [("Sleep apnoea", "3 definitions"), ("Insomnia", "2"), ("Constipation", "2"), ("ADHD", "2"), ("Epilepsy", "7")]
CAP = 8.0


def main():
    vs.apply()
    h = pd.read_csv(ROOT / "results" / "definition_effect" / "heterogeneity.tsv", sep="\t")
    cell = {(r.shared, r.family): r for r in h.itertuples()}
    cmap = vs.cmap_seq()
    norm = plt.Normalize(0, CAP)

    fig = plt.figure(figsize=(vs.SINGLE, 2.75))
    ax = fig.add_axes([0.30, 0.22, 0.68, 0.70])
    for j, t in enumerate(TRAITS):
        for i, (fam, _) in enumerate(SETS):
            r = cell[(t, fam)]
            lp = min(-np.log10(max(r.p, 1e-300)), CAP)
            size = 25 + 360 * r.rg_range
            ax.scatter(j, i, s=size, c=[cmap(norm(lp))], edgecolors=vs.INK if r.p < 0.05 else "none",
                       linewidths=0.8, zorder=3)
            if r.p < 0.05:
                ax.text(j, i + 0.47, f"Q {r.Q:.0f}, df {r.df}", ha="center", va="center", fontsize=7, color=vs.INK_2)
    for i, (fam, nd) in enumerate(SETS):
        ax.text(-0.6, i, f"{fam} ({nd.split()[0]})", ha="right", va="center", fontsize=7.5, color=vs.INK)
    for j, t in enumerate(TRAITS):
        ax.text(j, -0.75, TRAIT_SHORT[t], ha="center", va="center", fontsize=7.5, color=vs.INK_2)
        ax.scatter(j, -0.45, s=12, color=vs.TRAIT[t], clip_on=False)
    ax.set_xlim(-0.5, len(TRAITS) - 0.5)
    ax.set_ylim(len(SETS) - 0.3, -0.95)
    ax.axis("off")

    cax = fig.add_axes([0.08, 0.085, 0.36, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([0, 2, 4, 6, 8])
    cb.set_ticklabels(["0", "2", "4", "6", "8+"])
    cb.ax.tick_params(labelsize=7, colors=vs.MUTED, length=2)
    cb.outline.set_visible(False)
    cax.set_title("-log10 p of Q", fontsize=7, color=vs.INK_2, pad=2)
    kax = fig.add_axes([0.52, 0.045, 0.46, 0.09])
    kax.set_xlim(0, 3.8)
    kax.set_ylim(0, 1)
    kax.axis("off")
    kax.text(0.0, 0.5, "rg range", va="center", fontsize=7, color=vs.INK_2)
    for x, rr in ((1.05, 0.05), (2.0, 0.25), (3.05, 0.45)):
        kax.scatter(x, 0.5, s=25 + 360 * rr, color=vs.SEQ[3], edgecolors="none")
        kax.text(x + 0.30, 0.5, f"{rr:.2f}", va="center", fontsize=7, color=vs.INK_2)
    fig.text(0.02, 0.965, "Rows: condition (number of definitions). Ring: p < 0.05, labelled with Q and df.", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig13_heterogeneity")


if __name__ == "__main__":
    main()
