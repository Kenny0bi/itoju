"""Figure 8. SNP heritability of every definition, beside how many people it counts as cases.

A definition that casts a wider net gains cases, and with them precision, but it can also pull in
people whose condition is less genetic. Rows are grouped by condition, one row per definition.
Left: liability-scale SNP heritability with a 95% interval. Right: the number of cases on a log
scale, so a reader can see whether a definition with more people shows lower heritability.
A hollow mark at zero is an estimate below zero.

The liability conversion uses the sample prevalence as the population prevalence, a stated
approximation (FinnGen is a biobank, not a random population sample), so compare heights within a
condition rather than across conditions.

Source: results/ldsc/h2.tsv (standard two-step LDSC) and the FinnGen R12 manifest.
"""
import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs
from itoju_labels import GROUPS

ROOT = vs.ROOT
SETS = ["Epilepsy", "Sleep apnoea", "Insomnia", "Constipation", "ADHD", "Intellectual disability"]
Z = 1.96


def main():
    vs.apply()
    h2 = pd.read_csv(ROOT / "results" / "ldsc" / "h2.tsv", sep="\t").set_index("trait")
    man = pd.read_csv(ROOT / "data" / "raw" / "finngen_R12_manifest.tsv", sep="\t").set_index("phenocode")
    groups = dict(GROUPS)

    rows, headers, y = [], [], 0.0
    for name in SETS:
        headers.append((name, y))
        y += 0.95
        for e, lab in groups[name]:
            rows.append((e, lab, y))
            y += 1.0
        y += 0.4
    total = y

    fig = plt.figure(figsize=(vs.SINGLE, 0.16 * total + 0.85))
    a = fig.add_axes([0.43, 0.125, 0.29, 0.80])
    b = fig.add_axes([0.76, 0.125, 0.22, 0.80])
    for e, lab, yy in rows:
        r = h2.loc[e]
        cases = int(man.loc[e, "num_cases"])
        if r.h2 > 0:
            a.plot([max(r.h2 - Z * r.h2_se, 0), r.h2 + Z * r.h2_se], [yy, yy], color=vs.INK_2, lw=0.9, zorder=2)
            a.scatter(r.h2, yy, s=18, color=vs.INK, zorder=3)
        else:
            a.scatter(0, yy, s=18, facecolor="white", edgecolor=vs.INK, linewidths=0.8, zorder=3, clip_on=False)
        a.text(-0.022, yy, lab, ha="right", va="center", fontsize=7, color=vs.INK_2)
        b.plot([500, cases], [yy, yy], color=vs.SEQ[8], lw=2.2, solid_capstyle="butt")
        b.text(cases * 1.2, yy, f"{cases / 1000:.1f}k", va="center", fontsize=7, color=vs.INK_2)
    for name, yy in headers:
        a.text(-0.022, yy, name, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK)
        for ax, x0, x1 in ((a, 0, 0.33), (b, 500, 1e6)):
            ax.plot([x0, x1], [yy + 0.47, yy + 0.47], color=vs.GRID, lw=0.5)
    for ax in (a, b):
        ax.set_ylim(total, -0.8)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
    a.set_xlim(0, 0.33)
    a.set_xticks([0, 0.1, 0.2, 0.3])
    a.set_xlabel("liability h2, 95% interval")
    b.set_xscale("log")
    b.set_xlim(500, 1e6)
    b.set_xticks([1e3, 1e4, 1e5])
    b.set_xticklabels(["1k", "10k", "100k"])
    b.set_xlabel("cases, log scale")
    a.set_title("Heritability", loc="left", fontsize=8, color=vs.INK_2, pad=3)
    b.set_title("Cases", loc="left", fontsize=8, color=vs.INK_2, pad=3)
    fig.text(0.02, 0.006, "Hollow mark at zero: estimate below zero.", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig08_heritability")


if __name__ == "__main__":
    main()
