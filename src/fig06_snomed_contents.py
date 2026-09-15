"""Figure 6. What two definitions of one condition contain, and how much they share, counted in SNOMED CT.

Each definition's ICD-10 codes are translated to standard SNOMED CT concepts (OMOP "Maps to"). For each
pair of definitions, every concept either definition admits is one square, so set size and overlap are
read as area and position:
  ochre square   only in the first (narrower) definition
  split square   in both definitions
  blue square    only in the second definition
The ochre bracket above spans the first definition's concepts and the blue bracket below the second's,
each with its count. Squares sit in the order only-first, shared, only-second, so the shared part is
always where the two brackets meet.

Laxative purchases in the constipation definition have no SNOMED CT condition concept and are not drawn.
The published tables carry concept IDs only; the one annotated concept list is a short excerpt.
Sources: results/omop/endpoint_snomed_concepts.tsv (src/08b), data/ref/omop/target_concepts.tsv.
"""
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon, Rectangle

import viz_style as vs

ROOT = vs.ROOT
FIRST = "#c9761a"    # validated pair with SECOND (dataviz validator, all checks pass)
SECOND = "#2f6fc4"
PAIRS = [  # condition, first endpoint, its label, second endpoint, its label
    ("Sleep apnoea", "G6_SLEEPAPNO", "sleep apnoea", "SLEEP", "any sleep disorder"),
    ("Insomnia", "F5_INSOMNIA", "F51.0 or G47.0", "KRA_PSY_SLEEP_NONORG_EXMORE", "all of F51"),
    ("Constipation", "K11_CONSTIPATION", "K59.0 or laxatives", "K11_OTHFUNC", "all of K59"),
    ("ADHD", "F5_ADHD", "F90.0", "KRA_PSY_HYPERKIN_EXMORE", "all of F90"),
    ("Intellectual disability", "F5_MILDRET", "mild (F70)", "KRA_PSY_MENTALRET_EXMORE", "any (F7)"),
    ("Epilepsy", "FE", "focal", "G6_EPLEPSY", "any epilepsy"),
    ("Epilepsy", "FE", "focal", "GE", "generalized"),
]
NOTE_FOR = "K11_OTHFUNC"
NOTE_CONCEPTS = ["Irritable bowel syndrome", "Functional diarrhea", "Neurogenic bowel"]
X0 = 1.28          # inches: where the squares start
PITCH = 0.128      # inches per square
SIZE = 0.108
ROW = 0.60
NOTE = 0.24
TOP = 0.62


def main():
    vs.apply()
    sc = pd.read_csv(ROOT / "results" / "omop" / "endpoint_snomed_concepts.tsv", sep="\t", dtype=str)
    tgt = pd.read_csv(ROOT / "data" / "ref" / "omop" / "target_concepts.tsv", sep="\t", dtype=str,
                      keep_default_na=False, quoting=3).set_index("concept_id")
    concepts = lambda e: set(sc[sc.endpoint == e].snomed_concept_id)

    rows = []
    for cond, a, alab, b, blab in PAIRS:
        A, B = concepts(a), concepts(b)
        rows.append((cond, alab, blab, len(A - B), len(A & B), len(B - A), b, B - A))
        print(f"  {cond:24s} {alab:>20s} vs {blab:<20s} only first {len(A - B):2d}  shared {len(A & B):2d}  only second {len(B - A):2d}")
    widest = max(o + s + t for _, _, _, o, s, t, _, _ in rows)
    assert X0 + widest * PITCH < vs.SINGLE - 0.03, "squares would run off the column"

    fig_h = TOP + len(rows) * ROW + NOTE + 0.05
    fig = plt.figure(figsize=(vs.SINGLE, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, vs.SINGLE)
    ax.set_ylim(fig_h, 0)
    ax.axis("off")

    def square(x, y, kind):
        if kind == "both":
            ax.add_patch(Polygon([(x, y), (x + SIZE, y), (x, y + SIZE)], closed=True, facecolor=FIRST, edgecolor="none"))
            ax.add_patch(Polygon([(x + SIZE, y), (x + SIZE, y + SIZE), (x, y + SIZE)], closed=True, facecolor=SECOND,
                                 edgecolor="none"))
        else:
            ax.add_patch(Rectangle((x, y), SIZE, SIZE, facecolor=FIRST if kind == "first" else SECOND, edgecolor="none"))

    # key
    ax.text(0.05, 0.14, "Each square is one SNOMED CT concept a definition admits.", fontsize=7.5, color=vs.INK,
            va="center")
    kx = 0.05
    for kind, lab in (("first", "only in the first definition"), ("both", "in both"),
                      ("second", "only in the second")):
        square(kx, 0.33, kind)
        ax.text(kx + SIZE + 0.05, 0.33 + SIZE / 2, lab, fontsize=7, color=vs.INK_2, va="center")
        kx += SIZE + 0.1 + 0.052 * len(lab)

    y = TOP
    for cond, alab, blab, o, s, t, bkey, b_only in rows:
        sq_y = y + 0.22
        ax.text(X0 - 0.1, sq_y + SIZE / 2 - 0.035, cond, fontsize=7.6, color=vs.INK, ha="right", va="center")
        share = f"shares {s} of {o + s + t}" if s else "shares none"
        ax.text(X0 - 0.1, sq_y + SIZE / 2 + 0.09, share, fontsize=7, color=vs.INK_2, ha="right", va="center")
        kinds = ["first"] * o + ["both"] * s + ["second"] * t
        for k, kind in enumerate(kinds):
            square(X0 + k * PITCH, sq_y, kind)
        # first definition: bracket above; second: bracket below
        a0, a1 = X0, X0 + (o + s - 1) * PITCH + SIZE
        b0, b1 = X0 + o * PITCH, X0 + (o + s + t - 1) * PITCH + SIZE
        ax.plot([a0, a0, a1, a1], [sq_y - 0.02, sq_y - 0.06, sq_y - 0.06, sq_y - 0.02], color=FIRST, lw=0.9,
                solid_capstyle="butt")
        ax.text(a0, sq_y - 0.09, f"{alab}: {o + s}", fontsize=7, color=vs.INK, va="bottom")
        yb = sq_y + SIZE
        ax.plot([b0, b0, b1, b1], [yb + 0.02, yb + 0.06, yb + 0.06, yb + 0.02], color=SECOND, lw=0.9,
                solid_capstyle="butt")
        ax.text(b0, yb + 0.085, f"{blab}: {s + t}", fontsize=7, color=vs.INK, va="top")
        y += ROW
        if bkey == NOTE_FOR:
            names = {tgt.loc[c, "concept_name"] for c in b_only}
            assert all(n in names for n in NOTE_CONCEPTS), "note must name concepts the definition really adds"
            ax.text(0.05, y - 0.02, "all of K59 adds irritable bowel syndrome, functional diarrhoea,\nneurogenic bowel and more",
                    fontsize=7, color=vs.INK_2, va="top", linespacing=1.15)
            y += NOTE
    vs.save(fig, "fig06_snomed_contents")


if __name__ == "__main__":
    main()
