"""Figure 16. Autism in focus: how large a definition effect could this study see, and how large was it?

One row per alternative definition, compared with the first definition of its condition, with
iPSYCH-PGC autism as the shared trait. Each row is a caliper:
  stem          the observed shift in rg (this definition minus the first)
  dark sleeve   plus or minus 1.96 SE of the difference: shifts the test cannot tell from chance
  pale sleeve   the smallest shift the design would catch 4 times in 5 (alpha 0.05, power 0.80)
A stem that escapes the dark sleeve is colored. The two columns on the right give the observed shift
and the catchable one, so a well-powered null (small pale sleeve, stem inside) reads at a glance.

Source: results/definition_effect/pairwise_delta.tsv.
"""
import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs
from itoju_labels import GROUPS

ROOT = vs.ROOT
SETS = ["Epilepsy", "Sleep apnoea", "Insomnia", "Constipation", "ADHD"]
ASD = vs.TRAIT["ASD"]
INNER = "#d8d0bf"
OUTER = "#efebe2"
GRAY = "#aaa69c"
XCAP = 0.6


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw = pw[pw.shared == "ASD"]
    groups = dict(GROUPS)

    rows, headers, y = [], [], 0.0
    for name in SETS:
        defs = groups[name]
        headers.append((name, defs[0][1], y))
        y += 1.05
        for e, lab in defs[1:]:
            r = pw[((pw.def1 == defs[0][0]) & (pw.def2 == e)) | ((pw.def1 == e) & (pw.def2 == defs[0][0]))].iloc[0]
            shift = r.delta if r.def1 == e else -r.delta
            rows.append((lab, y, shift, 1.96 * r.se_delta, r.min_detectable_delta))
            y += 1.0
        y += 0.45
    total = y - 0.45

    fig = plt.figure(figsize=(vs.SINGLE, 0.14 * total + 0.88))
    ax = fig.add_axes([0.385, 0.19, 0.36, 0.73])
    for lab, yy, shift, hw, mdd in rows:
        ax.hlines(yy, -min(mdd, XCAP), min(mdd, XCAP), colors=OUTER, linewidth=6.5, zorder=1)
        ax.hlines(yy, -min(hw, XCAP), min(hw, XCAP), colors=INNER, linewidth=6.5, zorder=2)
        clear = abs(shift) > hw
        ax.hlines(yy, 0, shift, colors=ASD, linewidth=1.5 if clear else 1.0, zorder=4)
        ax.scatter(shift, yy, s=18, facecolor=ASD if clear else "white", edgecolor=ASD, linewidths=1.0, zorder=5)
        ax.text(-XCAP - 0.05, yy, lab, ha="right", va="center", fontsize=7, color=vs.INK_2)
        shown = "0.00" if abs(shift) < 0.005 else f"{shift:+.2f}"   # no "-0.00"
        ax.text(1.22, yy, shown, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=7, color=vs.INK if clear else vs.INK_2)
        ax.text(1.62, yy, f"{mdd:.2f}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7,
                color=vs.INK_2)
    name_x = (0.02 - 0.385) / 0.36          # figure left margin, in axes units
    for name, first, yy in headers:
        ax.text(name_x, yy, f"{name}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.6,
                fontweight="bold", color=vs.INK)
        # the header row carries no marks, so the reference definition can sit over the gauge
        ax.text(-XCAP - 0.02, yy, f"vs {first[0].lower() + first[1:] if not first[:2].isupper() else first}",
                ha="left", va="center", fontsize=7, color=vs.MUTED)
    # zero line drawn row by row, so it never runs through a header's text
    ax.vlines([0] * len(rows), [yy - 0.5 for _, yy, *_ in rows], [yy + 0.5 for _, yy, *_ in rows], colors=vs.INK_2,
              lw=0.5, zorder=3)
    ax.set_ylim(total, -1.1)
    ax.set_xlim(-XCAP - 0.02, XCAP + 0.02)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-0.5, 0, 0.5])
    ax.set_xticklabels(["-0.5", "0", "+0.5"])
    ax.set_xlabel("shift in rg with autism", fontsize=7.5)
    ax.text(1.22, -0.95, "shift", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7,
            color=vs.INK)
    ax.text(1.62, -0.95, "catchable", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7,
            color=vs.INK)
    # key
    kx = fig.add_axes([0.02, 0.005, 0.96, 0.075])
    kx.set_xlim(0, 100)
    kx.set_ylim(0, 1)
    kx.axis("off")
    kx.hlines(0.5, 0.5, 9.5, colors=OUTER, linewidth=6.5)
    kx.hlines(0.5, 3, 7, colors=INNER, linewidth=6.5)
    kx.text(11, 0.5, "Dark: shifts the test cannot tell from chance (1.96 SE).\n"
            "Pale: shift caught 4 times in 5. Filled dot: p < 0.05. Cut at 0.6.", va="center",
            fontsize=7, color=vs.INK_2, linespacing=1.3)
    vs.save(fig, "fig16_autism_focus")


if __name__ == "__main__":
    main()
