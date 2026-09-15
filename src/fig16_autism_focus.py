"""Figure 16. Autism in focus: does the way a co-occurring condition is defined change its genetic
overlap with autism?

One row per definition, grouped by condition. Left: rg with iPSYCH-PGC autism and its 95% interval
(standard two-step LDSC on full SNP sets). Right: for every definition after the first in its
condition, the observed shift from the first definition, drawn against the smallest shift this
design could detect for autism (gray bar; alpha 0.05, power 0.80, jackknife SE of the difference).
A dot inside its gray bar is a shift the data cannot tell from zero; a short gray bar with the dot
inside it is a well-powered null.

Sources: results/ldsc/rg.tsv, results/definition_effect/pairwise_delta.tsv.
"""
import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs
from itoju_labels import GROUPS

ROOT = vs.ROOT
SETS = ["Epilepsy", "Sleep apnoea", "Insomnia", "Constipation", "ADHD"]
Z = 1.96
ASD = vs.TRAIT["ASD"]


def main():
    vs.apply()
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    rg = rg[rg.trait1 == "ASD"].set_index("trait2")
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw = pw[pw.shared == "ASD"]
    groups = dict(GROUPS)

    rows, headers, y = [], [], 0.0
    for name in SETS:
        headers.append((name, y))
        y += 0.95
        defs = groups[name]
        for k, (e, lab) in enumerate(defs):
            rows.append((name, defs[0][0], e, lab, y, k == 0))
            y += 1.0
        y += 0.4
    total = y

    fig = plt.figure(figsize=(vs.SINGLE, 0.16 * total + 0.95))
    a = fig.add_axes([0.43, 0.155, 0.30, 0.765])
    b = fig.add_axes([0.77, 0.155, 0.21, 0.765])
    for name, ref, e, lab, yy, first in rows:
        r = rg.loc[e]
        a.plot([r.rg - Z * r.se, r.rg + Z * r.se], [yy, yy], color=ASD, lw=1.0, alpha=0.6)
        a.scatter(r.rg, yy, s=18, color=ASD, zorder=3)
        a.text(-0.72, yy, lab, ha="right", va="center", fontsize=7, color=vs.INK_2)
        if first:
            b.text(0, yy, "reference", ha="center", va="center", fontsize=7, color=vs.MUTED, zorder=5,
                   bbox=dict(boxstyle="square,pad=0.1", facecolor=vs.SURFACE, edgecolor="none"))
            continue
        row = pw[((pw.def1 == ref) & (pw.def2 == e)) | ((pw.def1 == e) & (pw.def2 == ref))].iloc[0]
        shift = row.delta if row.def1 == e else -row.delta   # this definition minus the reference
        m = row.min_detectable_delta
        b.fill_between([-m, m], yy - 0.3, yy + 0.3, color=vs.GRID, lw=0)
        b.scatter(shift, yy, s=16, color=ASD, zorder=3, edgecolors=vs.INK if row.p < 0.05 else "none", linewidths=0.7)
    for name, yy in headers:
        a.text(-0.72, yy, name, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK)
        for ax, x0, x1 in ((a, -0.7, 1.0), (b, -0.9, 0.9)):
            ax.plot([x0, x1], [yy + 0.47, yy + 0.47], color=vs.GRID, lw=0.5)
    for ax in (a, b):
        ax.set_ylim(total, -0.8)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.axvline(0, color=vs.MUTED, lw=0.6, zorder=0)
    a.set_xlim(-0.7, 1.0)
    a.set_xlabel("rg with autism\n(95% interval)")
    b.set_xlim(-0.9, 0.9)
    b.set_xticks([-0.5, 0, 0.5])
    b.set_xlabel("shift from\nfirst definition")
    a.set_title("Genetic correlation", loc="left", fontsize=8, color=vs.INK_2, pad=3)
    b.set_title("Definition shift", loc="left", fontsize=8, color=vs.INK_2, pad=3)
    fig.text(0.02, 0.012, "Gray bar: smallest shift detectable (80% power). Ringed dot: p < 0.05.", fontsize=7,
             color=vs.INK_2)
    vs.save(fig, "fig16_autism_focus")


if __name__ == "__main__":
    main()
