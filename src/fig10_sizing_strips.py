"""Figure 10. Sizing strips: the same condition name, measured through different definitions.

The form borrows from a clothing size chart: one label, several measurements. For each definition
set (sleep apnoea, insomnia, constipation, ADHD) and each psychiatric trait, every definition's rg
sits on one shared axis, joined by a thin line so the shift reads as a single movement. Behind
each row, a gray band spans the minimum detectable difference (alpha 0.05, power 0.80) around the
first definition: a second definition landing inside the band is a shift this design could not
tell from noise; landing outside it is a shift it would detect. Differences marked with a star
pass the jackknife test at p < 0.05.

Sources: results/definition_effect/rg_common.tsv and pairwise_delta.tsv (one-step LDSC on the
common SNP set, the estimates the difference test is built on).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs
from itoju_labels import ENDPOINT_LABEL, TRAIT_LABEL

ROOT = vs.ROOT
SETS = [
    ("Sleep apnoea", ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"]),
    ("Insomnia", ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"]),
    ("Constipation", ["K11_CONSTIPATION", "K11_OTHFUNC"]),
    ("ADHD", ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"]),
]
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
SHAPES = ["o", "D", "s"]


def main():
    vs.apply()
    rc = pd.read_csv(ROOT / "results" / "definition_effect" / "rg_common.tsv", sep="\t")
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    rg = {(r.trait1, r.trait2): (r.rg, r.se) for r in rc.itertuples()}

    fig, axes = plt.subplots(2, 2, figsize=(vs.DOUBLE, 4.6), sharex=True)
    fig.subplots_adjust(left=0.135, right=0.99, top=0.86, bottom=0.10, hspace=0.62, wspace=0.10)
    for k, (ax, (name, defs)) in enumerate(zip(axes.flat, SETS)):
        for i, t in enumerate(TRAITS):
            y = len(TRAITS) - 1 - i
            ref = defs[0]
            mdd = pw[(pw.shared == t) & (pw.def1 == ref) & (pw.def2 == defs[1])].min_detectable_delta.iloc[0]
            r0 = rg[(t, ref)][0]
            ax.fill_between([r0 - mdd, r0 + mdd], y - 0.28, y + 0.28, color=vs.GRID, lw=0, zorder=1)
            # each definition gets a small vertical offset inside its row, so identical estimates
            # (the finding in sleep apnoea and ADHD) stay visible instead of hiding one another
            offs = [0.17, 0.0, -0.17] if len(defs) == 3 else [0.11, -0.11]
            vals = [rg[(t, d)][0] for d in defs]
            ax.plot(vals, [y + o for o in offs], color=vs.TRAIT[t], lw=1.0, zorder=2, solid_capstyle="butt")
            for j, d in enumerate(defs):
                v, se = rg[(t, d)]
                yy = y + offs[j]
                ax.plot([v - se, v + se], [yy, yy], color=vs.TRAIT[t], lw=0.6, alpha=0.5, zorder=2)
                ax.scatter(v, yy, marker=SHAPES[j], s=30, facecolor=vs.TRAIT[t] if j == 0 else "white",
                           edgecolor=vs.TRAIT[t], linewidths=1.1, zorder=4)
                if j > 0:
                    row = pw[(pw.shared == t) & (((pw.def1 == ref) & (pw.def2 == d)) | ((pw.def1 == d) & (pw.def2 == ref)))]
                    if len(row) and row.p.iloc[0] < 0.05:
                        ax.text(v + se + 0.012, yy, "*", ha="left", va="center", fontsize=8.5, color=vs.INK)
            if k % 2 == 0:
                ax.text(-0.02, y, TRAIT_LABEL[t], transform=ax.get_yaxis_transform(), ha="right", va="center",
                        fontsize=7.5, color=vs.INK_2)
        ax.set_ylim(-0.6, len(TRAITS) - 0.3)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.axvline(0, color=vs.GRID, lw=0.6, zorder=0)
        ax.set_title(name, loc="left", fontsize=8.5, color=vs.INK, pad=16)
        # key for this panel's definitions
        for j, d in enumerate(defs):
            kx = 0.0 + j * (1.0 / len(defs))
            ax.scatter(kx + 0.015, 1.045, marker=SHAPES[j], s=22, transform=ax.transAxes, clip_on=False,
                       facecolor=vs.INK_2 if j == 0 else "white", edgecolor=vs.INK_2, linewidths=1.0)
            ax.text(kx + 0.04, 1.045, ENDPOINT_LABEL[d], transform=ax.transAxes, va="center", fontsize=7,
                    color=vs.INK_2)
    for ax in axes[1]:
        ax.set_xlabel("genetic correlation with the psychiatric trait, rg")
    axes[1, 0].set_xlim(-0.15, 0.95)
    fig.text(0.135, 0.965, "Gray band: smallest shift from the first definition this design detects (80% power).",
             fontsize=7, color=vs.INK_2)
    fig.text(0.135, 0.94, "* differs from the first definition, p < 0.05.   Thin bars: plus or minus one SE.",
             fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig10_sizing_strips")


if __name__ == "__main__":
    main()
