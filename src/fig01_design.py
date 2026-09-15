"""Figure 1. Study design, with the real numbers at every step.

Left to right: the data that came in, what was done to it, the analyses, and what each analysis
answers. Every count is from this study's own files and logs (FinnGen R12 manifest, GWAS file
headers, munge and LDSC logs, the OMOP mapping output).
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import viz_style as vs

ROOT = vs.ROOT


def box(ax, x, y, w, h, title, lines, face="white", edge=vs.MUTED, title_color=vs.INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.012",
                                facecolor=face, edgecolor=edge, linewidth=0.7))
    ax.text(x + 0.012, y + h - 0.022, title, fontsize=7.6, fontweight="bold", color=title_color, va="top")
    for i, line in enumerate(lines):
        ax.text(x + 0.012, y + h - 0.066 - i * 0.043, line, fontsize=7, color=vs.INK_2, va="top")


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=vs.MUTED, lw=0.7, shrinkA=0, shrinkB=0, mutation_scale=7))


def main():
    vs.apply()
    fig = plt.figure(figsize=(vs.DOUBLE, 3.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    col = [0.008, 0.258, 0.508, 0.758]
    w = 0.226
    heads = ["Data", "Preparation", "Analysis", "Question answered"]
    for x, h in zip(col, heads):
        ax.text(x, 0.975, h, fontsize=8.5, color=vs.INK, va="top")

    box(ax, col[0], 0.63, w, 0.30, "FinnGen R12 (Finland)",
        ["500,348 people",
         "28 endpoints: 6 definition",
         "   sets, 10 neighbours",
         "DF12 definition files",
         "LD matrix, 520,210 people"])
    box(ax, col[0], 0.29, w, 0.31, "Psychiatric GWAS (EUR)",
        ["Autism 18,381 / 27,969",
         "Schizophrenia 53,386 / 77,258",
         "Bipolar 41,917 / 371,549",
         "Depression 412,305 / 1.59M",
         "PTSD effective N 641,533",
         "  EHR 311,831; clinical 26,778"])
    # trait color dots beside the GWAS names (identity by mark, text stays ink)
    for i, t in enumerate(["ASD", "SCZ", "BIP", "MDD", "PTSD", "MDD_EHR"]):
        ax.scatter(col[0] + w - 0.012, 0.29 + 0.31 - 0.078 - i * 0.043, s=10, color=vs.TRAIT[t], zorder=3)
    box(ax, col[0], 0.05, w, 0.21, "Vocabularies, reference",
        ["OMOP: WHO ICD-10, ATC,",
         "   SNOMED CT, RxNorm",
         "1000 Genomes EUR LD scores"])

    box(ax, col[1], 0.63, w, 0.30, "Harmonise",
        ["keep HapMap3 SNPs",
         "1.17M SNPs per endpoint",
         "allele check vs 1000G EUR",
         "   r 0.991 to 0.998",
         "common set 953,474 SNPs"])
    box(ax, col[1], 0.29, w, 0.31, "Take definitions apart",
        ["expand regex code rules",
         "resolve INCLUDE chains",
         "sources: hospital, death,",
         "   primary care, drugs",
         "rules: exclusion, mode,",
         "   control sets"])
    box(ax, col[1], 0.05, w, 0.21, "Map to one vocabulary",
        ["153/153 ICD-10 to SNOMED",
         "99 concepts; ICD-8 unmapped",
         "phecodeX as second grouping"])

    box(ax, col[2], 0.63, w, 0.30, "LD score regression",
        ["autism h2 0.117 vs 0.118",
         "autism x ADHD 0.346 vs 0.360",
         "35 heritabilities",
         "217 genetic correlations",
         "intercepts absorb overlap"])
    box(ax, col[2], 0.29, w, 0.31, "Definition effect test",
        ["same 200 blocks per pair",
         "rg difference, jackknife SE",
         "135 comparisons, 5 sets",
         "heterogeneity Q, power",
         "31 EHR vs clinical checks"])
    box(ax, col[2], 0.05, w, 0.21, "Overlap and distance",
        ["ICD-10, phecode, SNOMED",
         "does distance predict drift?",
         "Finnish vs European LD"])

    box(ax, col[3], 0.63, w, 0.30, "Is the pipeline right?",
        ["reproduces published autism",
         "   results",
         "Finnish registry autism:",
         "   rg 0.73 to 0.98"], face=vs.SEQ[0])
    box(ax, col[3], 0.29, w, 0.31, "Does definition move rg?",
        ["which conditions and traits,",
         "by how much, and how small",
         "a shift the data could see"], face=vs.SEQ[0])
    box(ax, col[3], 0.05, w, 0.21, "Why, when it does",
        ["clinical scope (SNOMED)",
         "vs codes and data sources"], face=vs.SEQ[0])

    for yc in (0.78, 0.445, 0.155):
        for k in range(3):
            arrow(ax, col[k] + w + 0.004, yc, col[k + 1] - 0.004, yc)
    vs.save(fig, "fig01_design")


if __name__ == "__main__":
    main()
