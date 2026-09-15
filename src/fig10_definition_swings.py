"""Figure 10. Definition swings.

For each definition set and psychiatric trait, the genetic correlation under the first definition
is a filled dot. An arc swings from it to the estimate under each alternative definition (above the
track for the first alternative, below for the second); the arc's height grows with the distance
moved, so a definition that changes nothing draws an almost flat line. The shaded block under each
arc spans plus or minus 1.96 SE of that difference around the first estimate: an arc that lands
outside its block is a shift the jackknife test detects (colored); one inside is not (gray).

Sleep apnoea's shifts are small, so its panel uses a narrower rg axis (labelled).
The SE of each rg on its own is deliberately not drawn: two definitions share cases, so separate
intervals would overstate the noise in the difference. Source: rg_common.tsv, pairwise_delta.tsv.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath

import viz_style as vs
from itoju_labels import TRAIT_LABEL

ROOT = vs.ROOT

SLEEVE = "#e6e1d6"
GRAY = "#aaa69c"
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
SETS = [
    ("Sleep apnoea", ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"],
     ["hospital records", "+ primary care", "any sleep disorder"], (-0.03, 0.47), [0, 0.1, 0.2, 0.3, 0.4]),
    ("Insomnia", ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"], ["F51.0 or G47.0", "all of F51"],
     (-0.12, 0.92), [0, 0.2, 0.4, 0.6, 0.8]),
    ("Constipation", ["K11_CONSTIPATION", "K11_OTHFUNC"], ["K59.0 or laxatives", "all of K59"],
     (-0.12, 0.92), [0, 0.2, 0.4, 0.6, 0.8]),
    ("ADHD", ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"], ["F90.0", "all of F90"],
     (-0.12, 0.92), [0, 0.2, 0.4, 0.6, 0.8]),
]
H = 0.36
YLIM = (-0.55, len(TRAITS) - 0.35)


def arc(ax, x0, x1, y, side, color, lw, h, ls="-"):
    verts = [(x0, y), ((x0 + x1) / 2, y + side * 2 * h), (x1, y)]
    ax.add_patch(PathPatch(MPath(verts, [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]), fill=False,
                           edgecolor=color, lw=lw, ls=ls, capstyle="round", zorder=3))


def main():
    vs.apply()
    rc = pd.read_csv(ROOT / "results" / "definition_effect" / "rg_common.tsv", sep="\t")
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    rg = {(r.trait1, r.trait2): r.rg for r in rc.itertuples()}

    fig, axes = plt.subplots(2, 2, figsize=(vs.DOUBLE, 4.0))
    fig.subplots_adjust(left=0.115, right=0.985, top=0.775, bottom=0.11, hspace=0.74, wspace=0.08)
    fig.canvas.draw()
    for k, (ax, (name, defs, names, xlim, ticks)) in enumerate(zip(axes.flat, SETS)):
        bb = ax.get_window_extent()
        # data units per point in x and in y, so an arc's height can follow the distance it spans on paper
        x_per_pt = (xlim[1] - xlim[0]) / (bb.width * 72 / fig.dpi)
        y_per_pt = (YLIM[1] - YLIM[0]) / (bb.height * 72 / fig.dpi)
        for i, t in enumerate(TRAITS):
            y = len(TRAITS) - 1 - i
            r0 = rg[(t, defs[0])]
            ax.plot(xlim, [y, y], color=vs.GRID, lw=0.5, zorder=0)
            for j, d in enumerate(defs[1:]):
                side = 1 if j == 0 else -1
                row = pw[(pw.shared == t) & (pw.def1 == defs[0]) & (pw.def2 == d)].iloc[0]
                hw = 1.96 * row.se_delta
                r1 = rg[(t, d)]
                clear = abs(r1 - r0) > hw
                ax.add_patch(Rectangle((r0 - hw, y if side > 0 else y - H), 2 * hw, H, facecolor=SLEEVE,
                                       edgecolor="none", zorder=1))
                color = vs.TRAIT[t]
                h = min(H * 0.85, 0.45 * abs(r1 - r0) / x_per_pt * y_per_pt)
                arc(ax, r0, r1, y, side, color, 1.5 if clear else 0.9, h, "-" if clear else (0, (2.2, 1.6)))
                ax.scatter(r1, y, s=24, facecolor="white", edgecolor=color, linewidths=1.2, zorder=4)
            ax.scatter(r0, y, s=30, color=vs.TRAIT[t], linewidths=0, zorder=5)
            if k % 2 == 0:
                ax.text(-0.02, y, TRAIT_LABEL[t], transform=ax.get_yaxis_transform(), ha="right", va="center",
                        fontsize=7.5, color=vs.INK_2)
        ax.set_ylim(*YLIM)
        ax.set_xlim(*xlim)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{v:g}" for v in ticks])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.axvline(0, color=vs.MUTED, lw=0.5, zorder=0)
        # title placed in points above the subtitle, so it never collides whatever the panel height
        ax.annotate(name, (0, 1), xycoords="axes fraction", xytext=(0, 27 if len(names) == 3 else 17), textcoords="offset points",
                    fontsize=8.5, color=vs.INK, va="bottom")
        from_to = f"from {names[0]}  to  {names[1]}"
        if len(names) == 3:
            from_to = f"from {names[0]}  to  {names[1]} (above)\nor to {names[2]} (below)"
        ax.text(0.0, 1.02, from_to, transform=ax.transAxes, fontsize=7, color=vs.INK_2, va="bottom", linespacing=1.2)
        if k == 0:
            ax.text(1.0, 1.02, "narrower rg axis", transform=ax.transAxes, fontsize=7, color=vs.MUTED,
                    ha="right", va="bottom")
    for ax in axes[1]:
        ax.set_xlabel("genetic correlation with the psychiatric trait, rg")

    kax = fig.add_axes([0.115, 0.925, 0.87, 0.06])
    kax.set_xlim(0, 100)
    kax.set_ylim(0, 1)
    kax.axis("off")
    kax.add_patch(Rectangle((0.5, 0.2), 5, 0.55, facecolor=SLEEVE, edgecolor="none"))
    kax.scatter([3], [0.2], s=26, color=vs.TRAIT["BIP"], linewidths=0, zorder=3)
    arc(kax, 3, 9.5, 0.2, 1, vs.TRAIT["BIP"], 1.3, 0.4)
    kax.scatter([9.5], [0.2], s=22, facecolor="white", edgecolor=vs.TRAIT["BIP"], linewidths=1.2, zorder=3)
    kax.text(12, 0.62, "Filled dot: first definition.  Arc: the swing to another definition, higher when it moves further.",
             fontsize=7, color=vs.INK_2, va="center")
    kax.text(12, 0.12, "Shaded block: swings the test cannot tell from chance.  Solid arc: clears it (p < 0.05).  Dashed: does not.",
             fontsize=7, color=vs.INK_2, va="center")
    vs.save(fig, "fig10_definition_swings")


if __name__ == "__main__":
    main()
