"""Figure 12. Every definition shift, drawn inside its own noise.

One mark per comparison (166). The pale sleeve spans plus or minus 1.96 SE of the difference, the
band a shift must escape before the jackknife test calls it more than chance. The stem is the shift
itself: rg under the second definition minus rg under the first. A stem that escapes its sleeve is
drawn in the trait's color; one that stays inside is gray.

Families get different widths per comparison, so the five-comparison families stay readable next to
epilepsy's 105. Sleep apnoea's shifts are small and precise, so a magnified copy sits above them.
Source: results/definition_effect/pairwise_delta.tsv.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs
from itoju_labels import ENDPOINT_ORDER, TRAIT_LABEL

ROOT = vs.ROOT

SLEEVE = "#e6e1d6"
GRAY = "#b9b5ac"
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
YMAX = 0.5
FAMILIES = [  # family, title, width per comparison (pt)
    ("Sleep apnoea", "Sleep apnoea", 7.0),
    ("Insomnia", "Insomnia", 8.0),
    ("Constipation", "Constipation", 8.8),
    ("ADHD", "ADHD", 8.0),
    ("Epilepsy", "Epilepsy", 1.15),
    ("MDD definition (positive control)", "EHR vs clinical", 1.7),
]
PAIR_GAP = 0.8
AX_GAP_PT = 9.0
ZOOM_LIM = (-0.02, 0.078)


def ordered(sub, fam):
    """Rows in drawing order: pair by pair, traits in fixed order inside each pair."""
    if fam.startswith("MDD"):
        rank = {s: i for i, s in enumerate(TRAITS + ENDPOINT_ORDER)}
        return [list(sub.sort_values("shared", key=lambda s: s.map(rank)).itertuples())]
    groups = []
    for _, g in sub.groupby(["def1", "def2"], sort=False):
        g = g.set_index("shared").reindex([t for t in TRAITS if t in set(g.shared)]).reset_index()
        groups.append(list(g.itertuples()))
    return groups


def draw_marks(ax, X, R, fam, unit_pt, lim):
    hw = np.minimum([r.hw for r in R], lim)
    ax.vlines(X, -hw, hw, colors=SLEEVE, linewidth=max(0.8, 0.72 * unit_pt), capstyle="round", zorder=1)
    ax.axhline(0, color=vs.INK_2, lw=0.5, zorder=2)
    sh = np.array([r.shift for r in R])
    clear = np.array([r.clear for r in R])
    col = [vs.TRAIT["MDD" if fam.startswith("MDD") else r.shared] for r in R]
    widths = [max(0.5, 0.26 * unit_pt) if r.clear else max(0.35, 0.12 * unit_pt) for r in R]
    ax.vlines(X, 0, sh, colors=col, linewidth=widths, zorder=3)
    ax.scatter(X[clear], sh[clear], s=max(4, (0.55 * unit_pt) ** 2), c=[c for c, k in zip(col, clear) if k],
               linewidths=0, zorder=4)
    return clear


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw["shift"] = -pw.delta                     # delta = rg1 - rg2; draw second minus first
    pw["hw"] = 1.96 * pw.se_delta
    pw["clear"] = pw["shift"].abs() > pw.hw

    fig_w, fig_h = vs.DOUBLE, 3.55
    left_pt, right_pt = 44.0, 6.0
    layouts, widths = [], []
    for fam, title, unit in FAMILIES:
        groups = ordered(pw[pw.family == fam], fam)
        xs, x = [], 0.0
        for g in groups:
            for r in g:
                xs.append((x, r))
                x += 1
            x += PAIR_GAP
        span = x - PAIR_GAP
        layouts.append((fam, title, groups, xs, span))
        widths.append(span * unit + unit)
    avail = fig_w * 72 - left_pt - right_pt - AX_GAP_PT * (len(FAMILIES) - 1)
    scale = avail / sum(widths)

    fig = plt.figure(figsize=(fig_w, fig_h))
    bottom, height = 0.27, 0.49
    xpt = left_pt
    for k, ((fam, title, groups, xs, span), w) in enumerate(zip(layouts, widths)):
        wpt = w * scale
        ax = fig.add_axes([xpt / (fig_w * 72), bottom, wpt / (fig_w * 72), height])
        xpt += wpt + AX_GAP_PT
        X = np.array([p for p, _ in xs])
        R = [r for _, r in xs]
        unit_pt = wpt / (span + 1)
        clear = draw_marks(ax, X, R, fam, unit_pt, YMAX)
        xlim = (-0.8, span + 0.8)
        ax.set_xlim(*xlim)
        ax.set_ylim(-YMAX - 0.03, YMAX + 0.03)
        ax.set_xticks([])
        for s in ("bottom", "top", "right"):
            ax.spines[s].set_visible(False)
        if k == 0:
            ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
            ax.set_yticklabels(["-0.4", "-0.2", "0", "+0.2", "+0.4"])
            ax.set_ylabel("change in rg when the\ndefinition changes", fontsize=7.5)
        else:
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
        ax.text(0.0, 1.075, title, transform=ax.transAxes, fontsize=8, color=vs.INK, ha="left", va="bottom")
        ax.text(0.0, 1.005, f"{int(clear.sum())} of {len(R)} escape", transform=ax.transAxes, fontsize=7,
                color=vs.INK_2, ha="left", va="bottom")

        notes = []
        if fam == "Sleep apnoea":
            codes = {"G6_SLEEPAPNO": "A", "G6_SLEEPAPNO_INCLAVO": "B", "SLEEP": "C"}
            for g in groups:
                gx = np.mean([p for p, r in xs if (r.def1, r.def2) == (g[0].def1, g[0].def2)])
                ax.text(gx, -YMAX - 0.06, f"{codes[g[0].def1]} to {codes[g[0].def2]}", ha="center", va="top",
                        fontsize=7, color=vs.INK_2)
            notes = ["A hospital records", "B + primary care", "C any sleep disorder"]
            # magnified copy, set in from the y axis so its marks are never read against the main scale
            ins = ax.inset_axes([0.14, 0.62, 0.86, 0.34])
            draw_marks(ins, X, R, fam, unit_pt * 0.86, 1.0)
            ins.set_xlim(*xlim)
            ins.set_ylim(*ZOOM_LIM)
            ins.set_xticks([])
            ins.set_yticks([])
            ins.set_facecolor("#faf8f3")
            for s in ins.spines.values():
                s.set_color(vs.GRID)
                s.set_linewidth(0.5)
            ins.axhline(0.05, color=vs.GRID, lw=0.5, zorder=0)
            ins.text(xlim[0] + 0.2, 0.051, "+0.05", ha="left", va="bottom", fontsize=7, color=vs.MUTED)
            factor = (2 * (YMAX + 0.03) / height) / ((ZOOM_LIM[1] - ZOOM_LIM[0]) / (0.34 * height))
            ax.texts[1].set_text(f"{ax.texts[1].get_text()}, inset {factor:.0f}x")
        elif fam == "Insomnia":
            notes = ["F51.0 to", "all of F51"]
        elif fam == "Constipation":
            notes = ["K59.0 or", "laxatives to", "all of K59"]
        elif fam == "ADHD":
            notes = ["F90.0 to", "all of F90"]
        elif fam == "Epilepsy":
            sub = pw[pw.family == fam]
            kind = lambda d: "F" if d.startswith("FE") else ("G" if d.startswith("GE") else "A")
            k1, k2 = sub.def1.map(kind), sub.def2.map(kind)
            fg = (k1 != k2) & (k1 != "A") & (k2 != "A")
            same = k1 == k2
            anyv = ~fg & ~same
            notes = [f"focal vs generalized: {int(sub[fg].clear.sum())} of {int(fg.sum())}",
                     f"any vs a subtype: {int(sub[anyv].clear.sum())} of {int(anyv.sum())}",
                     f"strict or mode vs plain: {int(sub[same].clear.sum())} of {int(same.sum())}"]
        else:
            notes = ["depression,", "two definitions,", "31 comparisons"]
        y0 = -YMAX - (0.19 if fam == "Sleep apnoea" else 0.06)
        ax.text(xlim[0], y0, "\n".join(notes), ha="left", va="top", fontsize=7, color=vs.INK_2, linespacing=1.25)

    # key: traits on the first row, how to read a mark on the second
    kax = fig.add_axes([left_pt / (fig_w * 72), 0.865, 1 - (left_pt + right_pt) / (fig_w * 72), 0.125])
    kax.set_xlim(0, 100)
    kax.set_ylim(0, 1)
    kax.axis("off")
    x = 0.0
    for t in TRAITS:
        kax.scatter(x + 0.6, 0.8, s=22, color=vs.TRAIT[t], linewidths=0)
        kax.text(x + 1.6, 0.8, TRAIT_LABEL[t], va="center", fontsize=7.5, color=vs.INK_2)
        x += 3.0 + 1.25 * len(TRAIT_LABEL[t])
    kax.vlines([0.6, 2.4], [0.0, 0.08], [0.42, 0.34], colors=SLEEVE, linewidth=4.5, capstyle="round")
    kax.vlines([0.6], [0.21], [0.5], colors=vs.TRAIT["BIP"], linewidth=1.1)
    kax.scatter([0.6], [0.5], s=10, color=vs.TRAIT["BIP"], linewidths=0)
    kax.vlines([2.4], [0.21], [0.29], colors=vs.TRAIT["BIP"], linewidth=0.5)
    kax.text(4.2, 0.24, "Sleeve: shifts the test cannot tell from chance (1.96 SE of the difference); sleeves wider "
             "than 0.5 are cut at the edge.\nStem: the shift itself, in its trait's color; thick with a dot when it escapes its sleeve (p < 0.05).",
             va="center", fontsize=7, color=vs.INK_2, linespacing=1.3)
    vs.save(fig, "fig12_noise_skyline")


if __name__ == "__main__":
    main()
