"""Figure 1. Study design, drawn in the study's own visual language.

Four steps, left to right, each drawn with real numbers from one worked example (sleep apnoea and
depression), and each introducing a mark the later figures reuse:
  1  One condition, several borders: the three FinnGen R12 sleep apnoea definitions as nested discs,
     area proportional to cases.
  2  Each border meets five psychiatric GWAS in LD score regression.
  3  The genetic correlation swings when the border moves: depression with sleep apnoea under hospital
     records and under any sleep disorder (the arc of Fig. 10).
  4  Is the swing bigger than its noise? The shift against plus or minus 1.96 SE of the difference from
     the shared-block jackknife (the sleeve of Figs. 12, 15 and 16).
The strip underneath lists the checks that run beside the main design.

Sources: data/raw/finngen_R12_manifest.tsv, results/definition_effect/pairwise_delta.tsv, itoju_labels.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MPath

import viz_style as vs
from itoju_labels import ENDPOINT_ORDER, TRAIT_LABEL

ROOT = vs.ROOT
APNOEA = [("G6_SLEEPAPNO", "hospital records"), ("G6_SLEEPAPNO_INCLAVO", "+ primary care"),
          ("SLEEP", "any sleep disorder")]
TRAITS = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
SLEEVE = "#e6e1d6"
MDD = vs.TRAIT["MDD"]
DEF_RAMP = ["#a3d9d2", "#4fb4a8", "#12998c"]   # narrow to broad sleep apnoea definition
FIG_H = 3.05


def stage_axes(fig, k):
    x0 = 0.015 + k * 0.248
    ax = fig.add_axes([x0, 0.30, 0.215, 0.50])
    ax.axis("off")
    return ax


def head(fig, k, title, sub):
    x0 = 0.015 + k * 0.248
    fig.text(x0, 0.955, title, fontsize=8.5, color=vs.INK, va="top")
    fig.text(x0, 0.875, sub, fontsize=7, color=vs.INK_2, va="top")


def main():
    vs.apply()
    man = pd.read_csv(ROOT / "data" / "raw" / "finngen_R12_manifest.tsv", sep="\t").set_index("phenocode")
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    row = pw[(pw.shared == "MDD") & (pw.def1 == "G6_SLEEPAPNO") & (pw.def2 == "SLEEP")].iloc[0]
    n_endpoints = len(ENDPOINT_ORDER)
    fig = plt.figure(figsize=(vs.DOUBLE, FIG_H))

    # 1 the three borders as bars of cases on one scale; the dashed line carries the hospital-records length
    # down across the other two bars, so the extra cases read as the part past it
    ax = stage_axes(fig, 0)
    head(fig, 0, "1  One condition, 3 borders", "sleep apnoea in FinnGen R12, cases")
    cases = [int(man.loc[e, "num_cases"]) for e, _ in APNOEA]
    ax.set_xlim(0, max(cases) * 1.32)
    ax.set_ylim(3.0, -0.75)
    for i, ((e, lab), n) in enumerate(zip(APNOEA, cases)):
        ax.add_patch(FancyBboxPatch((0, i - 0.22), n, 0.44, boxstyle="round,pad=0,rounding_size=0.0",
                                    facecolor=DEF_RAMP[i], edgecolor="none"))
        ax.text(0, i - 0.3, lab, fontsize=7, color=vs.INK, va="bottom")
        ax.text(n + max(cases) * 0.02, i, f"{n:,}", fontsize=7, color=vs.INK, va="center", ha="left")
    ax.plot([cases[0], cases[0]], [0.7, 2.3], color=vs.INK_2, lw=0.7, ls=(0, (2, 2)))
    ax.text(cases[0], 2.55, "hospital-records length", fontsize=7, color=vs.INK_2, ha="center", va="center")

    # 2 three borders meet five GWAS
    ax = stage_axes(fig, 1)
    head(fig, 1, "2  Each meets 5 GWAS", f"LD score regression; {n_endpoints} endpoints in all")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.4, 4.4)
    for i, lab in enumerate(["hospital", "+ primary care", "any sleep"]):
        yl = 3.3 - i * 1.3
        for j, t in enumerate(TRAITS):
            ax.plot([0.44, 0.66], [yl, 4 - j], color=vs.GRID, lw=0.7, zorder=1)
        ax.add_patch(FancyBboxPatch((0.0, yl - 0.32), 0.44, 0.64, boxstyle="round,pad=0,rounding_size=0.08",
                                    facecolor=DEF_RAMP[i], edgecolor="none", lw=0.6, zorder=2))
        ax.text(0.22, yl, lab, ha="center", va="center", fontsize=7, color="white" if i == 2 else vs.INK, zorder=3)
    for j, t in enumerate(TRAITS):
        ax.scatter(0.68, 4 - j, s=26, color=vs.TRAIT[t], zorder=3)
        ax.text(0.73, 4 - j, TRAIT_LABEL[t], va="center", fontsize=7, color=vs.INK_2)

    # 3 the swing
    ax = stage_axes(fig, 2)
    head(fig, 2, "3  The correlation swings", "depression with sleep apnoea")
    r0, r1 = row.rg1, row.rg2
    lo, hi = 0.26, 0.51
    ax.set_xlim(lo, hi)
    ax.set_ylim(-1.3, 1.6)
    ax.plot([lo + 0.005, hi - 0.005], [0, 0], color=vs.GRID, lw=0.8)
    verts = [(r0, 0), ((r0 + r1) / 2, 1.5), (r1, 0)]
    ax.add_patch(PathPatch(MPath(verts, [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]), fill=False, edgecolor=MDD, lw=1.6))
    ax.scatter(r0, 0, s=40, color=MDD, zorder=3)
    ax.scatter(r1, 0, s=34, facecolor="white", edgecolor=MDD, linewidths=1.4, zorder=3)
    ax.text(r0 + 0.004, -0.35, f"hospital records\nrg {r0:.2f}", ha="right", va="top", fontsize=7, color=vs.INK_2,
            linespacing=1.2)
    ax.text(r1 - 0.004, -0.35, f"any sleep disorder\nrg {r1:.2f}", ha="left", va="top", fontsize=7, color=vs.INK_2,
            linespacing=1.2)

    # 4 the sleeve
    ax = stage_axes(fig, 3)
    head(fig, 3, "4  Bigger than its noise?", "shared-block jackknife difference test")
    shift, hw = r1 - r0, 1.96 * row.se_delta
    ax.set_xlim(-0.03, 0.095)
    ax.set_ylim(-1.3, 1.6)
    ax.plot([-hw, hw], [0.3, 0.3], color=SLEEVE, lw=11, solid_capstyle="butt")
    ax.plot([0, 0], [-0.15, 0.75], color=vs.INK_2, lw=0.6)
    ax.plot([0, shift], [0.3, 0.3], color=MDD, lw=1.8)
    ax.scatter(shift, 0.3, s=34, color=MDD, zorder=3)
    ax.text(0, 1.0, f"noise: 1.96 SE = {hw:.3f}", ha="left", va="bottom", fontsize=7, color=vs.INK_2)
    mant, expo = f"{row.p:.0e}".split("e")
    ax.text(shift, -0.35, f"shift {shift:+.3f}\n" + rf"p = {mant} $\times$ 10$^{{{int(expo)}}}$", ha="right", va="top",
            fontsize=7, color=vs.INK_2, linespacing=1.2)

    for k in range(3):
        x = 0.015 + k * 0.248 + 0.222
        fig.add_artist(plt.Line2D([x, x + 0.016], [0.55, 0.55], color=vs.MUTED, lw=0.8, transform=fig.transFigure))
        fig.add_artist(plt.Line2D([x + 0.009, x + 0.016, x + 0.009], [0.575, 0.55, 0.525], color=vs.MUTED, lw=0.8,
                                  transform=fig.transFigure))

    fig.add_artist(plt.Line2D([0.015, 0.985], [0.25, 0.25], color=vs.GRID, lw=0.6, transform=fig.transFigure))
    fig.text(0.015, 0.21, "Checks beside the design", fontsize=7.5, color=vs.INK, va="center")
    checks = ["published autism results reproduced", "positive control: EHR vs clinical depression",
              "codes mapped to SNOMED CT", "Finnish vs European LD"]
    # two rows of two, so every check reads in full
    for k, c in enumerate(checks):
        x, yk = (0.02, 0.50)[k % 2], (0.155, 0.11)[k // 2]
        fig.add_artist(plt.Line2D([x], [yk], marker="o", markersize=3.5, color=vs.INK_2, ls="",
                                  transform=fig.transFigure))
        fig.text(x + 0.009, yk, c, fontsize=7, color=vs.INK_2, va="center")
    fig.text(0.015, 0.045, f"Data: FinnGen R12 summary statistics and endpoint definitions ({n_endpoints} endpoints); "
             "GWAS of autism, schizophrenia,\nbipolar disorder, depression (with EHR and clinical subsets) and PTSD; "
             "OMOP standard vocabularies (Athena).", fontsize=7, color=vs.INK_2, va="center", linespacing=1.3)
    vs.save(fig, "fig01_design")


if __name__ == "__main__":
    main()
