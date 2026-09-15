"""Figure 12. Detectable or not: every definition comparison placed against its own power.

x is the minimum difference in rg the comparison could detect (alpha 0.05, power 0.80), from the
jackknife SE of the difference. y is the difference actually observed. The line |delta| =
1.96 SE (drawn as y = 0.70 x, since the minimum detectable difference is 2.80 SE) separates
differences that pass p < 0.05 (above) from those that do not (below). A point far to the left
and below the line is a well-powered null; far to the right is a comparison too weak to say
anything. The idea comes from evidence-strength charts that place each claim against a threshold
rather than printing a p-value beside it.

Panel a: the 135 comparisons between definitions of the same medical condition.
Panel b: the 31 comparisons between EHR and clinical depression (positive control).
Source: results/definition_effect/pairwise_delta.tsv
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
RATIO = norm.ppf(0.975) / (norm.ppf(0.975) + norm.ppf(0.80))


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    within = pw[pw.family != "MDD definition (positive control)"]
    control = pw[pw.family == "MDD definition (positive control)"]

    fig, (a, b) = plt.subplots(1, 2, figsize=(vs.DOUBLE, 3.45), gridspec_kw={"width_ratios": [1.35, 1]})
    fig.subplots_adjust(left=0.08, right=0.99, top=0.80, bottom=0.20, wspace=0.28)
    lim = 0.62
    for ax in (a, b):
        xs = np.linspace(0, 1.25, 50)
        ax.fill_between(xs, 0, RATIO * xs, color=vs.GRID, lw=0, zorder=0)
        ax.plot(xs, RATIO * xs, color=vs.MUTED, lw=0.8, zorder=1)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, 0.48)
        ax.set_xlabel("smallest detectable difference in rg")
    a.set_ylabel("observed |difference| in rg")

    clipped = 0
    for r in within.itertuples():
        x = min(r.min_detectable_delta, lim - 0.008)
        clipped += r.min_detectable_delta > lim
        a.scatter(x, abs(r.delta), marker=FAMILY_SHAPE[r.family], s=18 if r.family == "Epilepsy" else 26,
                  facecolor=vs.TRAIT[r.shared], edgecolor="white", linewidths=0.5, zorder=3,
                  alpha=0.9 if r.family != "Epilepsy" else 0.55)
    a.text(lim - 0.005, 0.475, f"{clipped} weaker comparisons (generalized\nepilepsy) drawn at the right edge", ha="right",
           va="top", fontsize=7, color=vs.INK_2)
    # name the clearest findings
    notes = [("MDD", "G6_SLEEPAPNO", "SLEEP", "sleep apnoea vs any sleep\ndisorder (all traits)", (0.012, 0.445)),
             ("PTSD", "K11_CONSTIPATION", "K11_OTHFUNC", "constipation vs all functional\nbowel disorders (PTSD)", (0.05, 0.365)),
             ("MDD", "F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE", "insomnia definitions\n(MDD, BIP)", (0.05, 0.285))]
    # plain text plus a separately drawn leader line: an annotation's extent includes its arrow,
    # which the layout check (rightly) reads as covering the neighbouring labels
    for t, d1, d2, txt, (tx, ty) in notes:
        r = within[(within.shared == t) & (within.def1 == d1) & (within.def2 == d2)].iloc[0]
        px, py = r.min_detectable_delta, abs(r.delta)
        a.text(tx, ty, txt, fontsize=7, color=vs.INK_2, va="center", ha="left")
        # the leader ends just below the label's first characters, so it never crosses the text
        a.plot([px - 0.004, tx + 0.004], [py + 0.006, ty - 0.05], color=vs.MUTED, lw=0.5, zorder=2)
    a.set_title("a  Definitions of the same condition (135 comparisons)", loc="left", fontsize=8.5, pad=6)

    for r in control.itertuples():
        x = min(r.min_detectable_delta, lim - 0.008)
        kind = "s" if r.shared in TRAITS else "o"
        col = vs.TRAIT[r.shared] if r.shared in TRAITS else vs.MUTED
        b.scatter(x, abs(r.delta), marker=kind, s=26, facecolor=col, edgecolor="white", linewidths=0.5, zorder=3)
    top = control.sort_values("p").iloc[0]
    b.text(top.min_detectable_delta + 0.06, abs(top.delta) + 0.06, "PTSD", fontsize=7, color=vs.INK_2, va="center")
    b.plot([top.min_detectable_delta + 0.006, top.min_detectable_delta + 0.054], [abs(top.delta), abs(top.delta) + 0.06],
           color=vs.MUTED, lw=0.5, zorder=2)
    b.set_title("b  EHR vs clinical depression (31 comparisons)", loc="left", fontsize=8.5, pad=6)
    n_edge = int((control.min_detectable_delta > lim).sum())
    # keys in two rows; spacing follows label length so no symbol lands on a neighbouring word
    fig.subplots_adjust(top=0.80)
    def key_row(items, y):
        x = 0.08
        for lab, style in items:
            fig.add_artist(plt.Line2D([x], [y], ls="", transform=fig.transFigure, **style))
            fig.text(x + 0.009, y, lab, fontsize=7, color=vs.INK_2, va="center")
            x += 0.022 + 0.0105 * len(lab)
    key_row([(TRAIT_LABEL[t], dict(marker="o", markersize=4.5, markerfacecolor=vs.TRAIT[t], markeredgecolor=vs.TRAIT[t]))
             for t in TRAITS], 0.965)
    shapes = [(fam if fam != "Epilepsy" else "Epilepsy (drawn lighter, 105 pairs)",
               dict(marker=m, markersize=4.5, markerfacecolor="white", markeredgecolor=vs.INK_2)) for fam, m in FAMILY_SHAPE.items()]
    key_row(shapes, 0.915)
    fig.text(0.08, 0.012, "Line: |difference| = 1.96 SE. Points above it differ at p < 0.05; the gray region is not significant.",
             fontsize=7, color=vs.INK_2, ha="left", va="bottom")
    b.text(0.02, 0.475, f"gray circles: FinnGen medical endpoints\ncolored squares: psychiatric traits\n{n_edge} comparisons drawn at the right edge", fontsize=7, color=vs.INK_2, va="top")
    vs.save(fig, "fig12_detectability")


if __name__ == "__main__":
    main()
