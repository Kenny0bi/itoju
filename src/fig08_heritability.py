"""Figure 8. Does a wider net lower heritability? SNP heritability against the number of cases.

A definition that casts a wider net gains cases, and with them precision, but it can also pull in
people whose condition is less genetic. One panel per condition; each definition is placed at its
number of cases (log scale) and its liability-scale SNP heritability. The dot is the estimate, the
thick bar its 50% interval and the thin line its 95% interval. All panels share the heritability
axis. A hollow mark on the floor is an estimate below zero.

The liability conversion uses the sample prevalence as the population prevalence, a stated
approximation (FinnGen is a biobank, not a random population sample), so compare heights within a
condition rather than across conditions.

Source: results/ldsc/h2.tsv (standard two-step LDSC) and the FinnGen R12 manifest.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs

ROOT = vs.ROOT
TEAL = "#12998c"
TEAL_DARK = "#0b6b61"
PANELS = [
    ("Epilepsy", [("G6_EPLEPSY", "any"), ("FE", "focal"), ("FE_STRICT", "focal, strict"), ("FE_MODE", "focal, mode"),
                  ("GE", "generalized"), ("GE_STRICT", "generalized, strict"), ("GE_MODE", "generalized, mode")]),
    ("Sleep apnoea", [("G6_SLEEPAPNO", "hospital records"), ("G6_SLEEPAPNO_INCLAVO", "+ primary care"),
                      ("SLEEP", "any sleep disorder")]),
    ("Insomnia", [("F5_INSOMNIA", "F51.0 or G47.0"), ("KRA_PSY_SLEEP_NONORG_EXMORE", "all of F51")]),
    ("Constipation", [("K11_CONSTIPATION", "K59.0 or laxatives"), ("K11_OTHFUNC", "all of K59")]),
    ("ADHD", [("F5_ADHD", "F90.0"), ("KRA_PSY_HYPERKIN_EXMORE", "all of F90")]),
    ("Intellectual disability", [("F5_MILDRET", "mild (F70)"), ("KRA_PSY_MENTALRET_EXMORE", "any (F7)")]),
]
YMAX = 0.30
Z50, Z95 = 0.674, 1.96
FIG_H = 4.4
W, H = 0.295, 0.33
TICKS = [(1e3, "1k"), (2e3, "2k"), (5e3, "5k"), (1e4, "10k"), (2e4, "20k"), (5e4, "50k"), (1e5, "100k")]


def main():
    vs.apply()
    h2 = pd.read_csv(ROOT / "results" / "ldsc" / "h2.tsv", sep="\t").set_index("trait")
    man = pd.read_csv(ROOT / "data" / "raw" / "finngen_R12_manifest.tsv", sep="\t").set_index("phenocode")
    fig = plt.figure(figsize=(vs.DOUBLE, FIG_H))
    for p, (name, defs) in enumerate(PANELS):
        col, row = p % 3, p // 3
        ax = fig.add_axes([0.075 + col * 0.318, 0.575 - row * 0.455, W, H])
        cases = np.array([man.loc[e, "num_cases"] for e, _ in defs], float)
        lc = np.log10(cases)
        mid = (lc.min() + lc.max()) / 2
        half = max(0.4, (lc.max() - lc.min()) / 2 + 0.15)
        lo, hi = mid - half, mid + half
        lab_w = (hi - lo) * 0.95
        ax.set_xlim(lo, hi + lab_w)
        ax.set_ylim(0, YMAX)
        pts = []
        for (e, lab), x in zip(defs, lc):
            r = h2.loc[e]
            if r.h2 > 0:
                ax.plot([x, x], [max(r.h2 - Z95 * r.h2_se, 0), min(r.h2 + Z95 * r.h2_se, YMAX)], color=TEAL,
                        lw=0.9, zorder=2)
                ax.plot([x, x], [max(r.h2 - Z50 * r.h2_se, 0), r.h2 + Z50 * r.h2_se], color=TEAL_DARK, lw=2.8,
                        solid_capstyle="butt", zorder=3)
                ax.scatter(x, r.h2, s=24, facecolor=TEAL, edgecolor="white", linewidths=0.8, zorder=4)
                pts.append((r.h2, x, lab))
            else:
                ax.scatter(x, 0.0, s=18, facecolor="white", edgecolor=vs.MUTED, linewidths=1.0, zorder=4, clip_on=False)
                pts.append((0.0, x, f"{lab}: below 0"))
        # labels stacked in the right column, nudged apart only as far as needed
        gap = 8.5 / (H * FIG_H * 72) * YMAX
        pts.sort()
        ys = []
        for yv, _, _ in pts:
            ys.append(max(yv, ys[-1] + gap) if ys else max(yv, gap * 0.4))
        over = ys[-1] - (YMAX - gap * 0.5)
        if over > 0:
            ys = [v - over for v in ys]
        x_text = hi + lab_w * 0.10
        for (yv, x, lab), yl in zip(pts, ys):
            ax.plot([x + 0.03 * (hi - lo), x_text - 0.03 * (hi - lo)], [yv, yl], color=vs.GRID, lw=0.6, zorder=1)
            ax.text(x_text, yl, lab, va="center", fontsize=7, color=vs.INK_2)
        ticks = [(np.log10(v), s) for v, s in TICKS if lo <= np.log10(v) <= hi]
        ax.set_xticks([t for t, _ in ticks])
        ax.set_xticklabels([s for _, s in ticks])
        ax.spines["bottom"].set_bounds(lo, hi)
        ax.set_yticks([0, 0.1, 0.2, 0.3])
        if col == 0:
            ax.set_ylabel("liability SNP heritability", fontsize=7.5)
        else:
            ax.set_yticklabels([])
        if row == 1:
            ax.set_xlabel("cases (log scale)", fontsize=7.5, loc="left")
        ax.text(0.0, 1.04, name, transform=ax.transAxes, fontsize=8.5, color=vs.INK, va="bottom")
    fig.text(0.075, 0.975, "Teal dot: estimate.  Thick bar: 50% interval.  Thin line: 95% interval.  "
             "Hollow gray dot on the floor: estimate below zero.", fontsize=7, color=vs.INK_2, va="top")
    vs.save(fig, "fig08_heritability")


if __name__ == "__main__":
    main()
