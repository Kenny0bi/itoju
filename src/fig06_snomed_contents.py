"""Figure 6. What a definition contains, read in SNOMED CT.

Codes hide scope. Translating each definition's ICD-10 codes to standard SNOMED CT concepts
(OMOP "Maps to") shows the clinical ideas a definition actually lets in. Each row is one SNOMED
concept; each column one definition; a filled dot means the definition includes it.
  a  Insomnia: F51.0 or G47.0 lets in 2 concepts; all of F51 adds parasomnias and hypersomnia.
  b  Constipation: the constipation code (or a laxative purchase) is 1 concept; all of K59, the
     FinnGen "other functional intestinal disorders" endpoint, adds irritable bowel syndrome,
     functional diarrhoea, neurogenic bowel and post-surgical conditions.
  c  ADHD: F90.0 is 2 concepts; F90 adds hyperkinetic conduct disorder and hyperkinetic
     syndrome with developmental delay.
Laxative purchases have no condition concept, so they appear only as a note in panel b.

Sources: results/omop/endpoint_snomed_concepts.tsv (src/08b) and data/ref/omop/target_concepts.tsv.
"""
import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs
from itoju_labels import ENDPOINT_LABEL

ROOT = vs.ROOT
PANELS = [
    ("a  Insomnia", ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"]),
    ("b  Constipation", ["K11_CONSTIPATION", "K11_OTHFUNC"]),
    ("c  ADHD", ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"]),
]
SHORT_COL = {"F5_INSOMNIA": "F51.0,\nG47.0", "KRA_PSY_SLEEP_NONORG_EXMORE": "F51",
             "K11_CONSTIPATION": "K59.0 or\nlaxatives", "K11_OTHFUNC": "K59",
             "F5_ADHD": "F90.0", "KRA_PSY_HYPERKIN_EXMORE": "F90"}


def main():
    vs.apply()
    sc = pd.read_csv(ROOT / "results" / "omop" / "endpoint_snomed_concepts.tsv", sep="\t", dtype=str)
    tgt = pd.read_csv(ROOT / "data" / "ref" / "omop" / "target_concepts.tsv", sep="\t", dtype=str,
                      keep_default_na=False, quoting=3).set_index("concept_id")

    heights = []
    tables = []
    for title, defs in PANELS:
        sub = sc[sc.endpoint.isin(defs)]
        members = sub.groupby("snomed_concept_id").endpoint.apply(set)
        rows = sorted(members.index, key=lambda c: (-len(members[c]), tgt.loc[c, "concept_name"].lower()))
        tables.append((title, defs, rows, members))
        heights.append(len(rows) + 2.6 + (2.4 if title.startswith("b") else 0))

    fig = plt.figure(figsize=(vs.SINGLE, 0.155 * sum(heights) + 0.3))
    top = 0.985
    for (title, defs, rows, members), h in zip(tables, heights):
        frac = h / sum(heights) * 0.96
        ax = fig.add_axes([0.02, top - frac, 0.96, frac])
        top -= frac
        ax.set_xlim(0, 1)
        ax.set_ylim(h - 0.3, -0.5)
        ax.axis("off")
        ax.text(0.0, 0.2, title, fontsize=8.5, color=vs.INK, va="center")
        xs = [0.78, 0.92]
        for x, d in zip(xs, defs):
            ax.text(x, 1.0, SHORT_COL[d], ha="center", va="center", fontsize=7, color=vs.INK_2)
        for i, c in enumerate(rows):
            y = i + 2.2
            name = tgt.loc[c, "concept_name"]
            if len(name) > 52:
                name = name.replace("Attention deficit hyperactivity disorder", "ADHD")
            both = len(members[c]) == 2
            ax.text(0.0, y, name, fontsize=7, color=vs.INK if both else vs.INK_2, va="center")
            ax.plot([0.0, 0.95], [y + 0.5, y + 0.5], color=vs.GRID, lw=0.3)
            for x, d in zip(xs, defs):
                inside = d in members[c]
                ax.scatter(x, y, s=26 if inside else 8, facecolor=vs.INK if inside else "white",
                           edgecolor=vs.INK if inside else vs.GRID, linewidths=0.6, zorder=3)
        if title.startswith("b"):
            ax.text(0.0, len(rows) + 2.6, "Laxative purchases (ATC A06A, 69 substances)\nhave no SNOMED condition concept.",
                    fontsize=7, color=vs.MUTED, va="top")
    vs.save(fig, "fig06_snomed_contents")


if __name__ == "__main__":
    main()
