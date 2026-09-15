"""Figure: what two definitions of one condition let in, counted in SNOMED CT concepts.

Each definition's ICD-10 codes are translated to standard SNOMED CT concepts through the OMOP "Maps to"
relationship. For each pair of definitions, every concept either one admits is one square, in the order
only-first, shared, only-second, so set size and overlap read as length and position. The count of
shared concepts is printed at the right, and the concepts the wider constipation definition adds are
named (a short excerpt; the published tables carry concept IDs only).

Sources: results/omop/endpoint_snomed_concepts.tsv (src/08b), data/ref/omop/target_concepts.tsv.
"""
import pandas as pd

import itoju_svg as sv

ROOT = sv.ROOT
PAIRS = [  # condition, first endpoint, its label, second endpoint, its label
    ("Sleep apnoea", "G6_SLEEPAPNO", "sleep apnoea", "SLEEP", "any sleep disorder"),
    ("Insomnia", "F5_INSOMNIA", "F51.0 or G47.0", "KRA_PSY_SLEEP_NONORG_EXMORE", "all of F51"),
    ("Constipation", "K11_CONSTIPATION", "K59.0 or laxatives", "K11_OTHFUNC", "all of K59"),
    ("ADHD", "F5_ADHD", "F90.0", "KRA_PSY_HYPERKIN_EXMORE", "all of F90"),
    ("Intellectual disability", "F5_MILDRET", "mild, F70", "KRA_PSY_MENTALRET_EXMORE", "any, F7"),
    ("Epilepsy", "FE", "focal", "G6_EPLEPSY", "any epilepsy"),
    ("Epilepsy", "FE", "focal", "GE", "generalized"),
]
NOTE_FOR = "K11_OTHFUNC"
NOTE_NAMES = ["Irritable bowel syndrome", "Functional diarrhea", "Neurogenic bowel"]
SQ, PITCH = 9.0, 11.5


def main():
    sc = pd.read_csv(ROOT / "results" / "omop" / "endpoint_snomed_concepts.tsv", sep="\t", dtype=str)
    tgt = pd.read_csv(ROOT / "data" / "ref" / "omop" / "target_concepts.tsv", sep="\t", dtype=str,
                      keep_default_na=False, quoting=3).set_index("concept_id")
    concepts = lambda e: set(sc[sc.endpoint == e].snomed_concept_id)
    rows = []
    for cond, a, alab, b, blab in PAIRS:
        A, B = concepts(a), concepts(b)
        rows.append(dict(cond=cond, alab=alab, blab=blab, only_a=len(A - B), both=len(A & B), only_b=len(B - A),
                         b=b, added=sorted(tgt.loc[c, "concept_name"] for c in B - A)))
        print(f"  {cond:24s} {alab:>18s} vs {blab:<18s} only first {len(A - B):2d} shared {len(A & B):2d} only second {len(B - A):2d}")

    W = sv.DOUBLE
    left, x_sq, x_share = 14.0, 176.0, 400.0
    row_h, top = 34.0, 104.0
    note_h = 12.0
    H = top + len(rows) * row_h + note_h + 78
    f = sv.Figure(W, H)
    f.header("itoju  /  what a definition lets in",
             "The same condition name can let in a handful of diagnoses or dozens",
             "Each square is one standard SNOMED CT concept a definition's ICD-10 codes map to (OMOP, Athena)")
    f.key_row(left, 70, [("square", sv.FIRST, "only in the first definition"),
                         ("square", sv.SECOND, "only in the second")])
    kx = left + 11 + sv.text_width("only in the first definition") + 10 + 11 + sv.text_width("only in the second") + 10
    f.polygon([(kx, 64), (kx + 6, 64), (kx, 70)], sv.FIRST)
    f.polygon([(kx + 6, 64), (kx + 6, 70), (kx, 70)], sv.SECOND)
    f.text(kx + 11, 70, "in both", 7.0, sv.INK_2)
    f.text(x_share, 92, "shared", 7.0, sv.DIM)

    y = top
    for r in rows:
        f.text(left, y + 4, r["cond"], 7.6, sv.INK)
        f.text(left, y + 14, f"{r['alab']}  vs  {r['blab']}", 7.0, sv.INK_2)
        kinds = ["a"] * r["only_a"] + ["ab"] * r["both"] + ["b"] * r["only_b"]
        assert x_sq + len(kinds) * PITCH < x_share - 8, "squares run into the shared column"
        for k, kind in enumerate(kinds):
            x0, y0 = x_sq + k * PITCH, y - 3
            if kind == "ab":
                f.polygon([(x0, y0), (x0 + SQ, y0), (x0, y0 + SQ)], sv.FIRST, mark=True)
                f.polygon([(x0 + SQ, y0), (x0 + SQ, y0 + SQ), (x0, y0 + SQ)], sv.SECOND, mark=True)
            else:
                f.rect(x0, y0, SQ, SQ, sv.FIRST if kind == "a" else sv.SECOND, rx=1, mark=True)
        n_a, n_b, total = r["only_a"] + r["both"], r["only_b"] + r["both"], len(kinds)
        f.text(x_sq, y + 17, f"{r['alab']} {n_a}, {r['blab']} {n_b}", 7.0, sv.DIM)
        f.text(x_share, y + 5, f"{r['both']} of {total}", 9.0, sv.INK, family=sv.SERIF)
        f.text(x_share + 34, y + 5, "none" if r["both"] == 0 else f"{100 * r['both'] / total:.0f}%", 7.0, sv.DIM)
        y += row_h
        if r["b"] == NOTE_FOR:
            assert all(n in r["added"] for n in NOTE_NAMES), r["added"]
            f.text(x_sq, y - 7, "adds irritable bowel syndrome, functional diarrhoea, neurogenic bowel and more", 7.0,
                   sv.INK_2)
            y += note_h

    # the reading, computed from the rows
    by = {(r["cond"], r["blab"]): r for r in rows}
    apn, con, epi = by[("Sleep apnoea", "any sleep disorder")], by[("Constipation", "all of K59")], by[("Epilepsy", "generalized")]
    lines = [
        f"Widening sleep apnoea to any sleep disorder adds {apn['only_b']} concepts to the {apn['both']} they share.",
        f"All of K59 adds {con['only_b']} to constipation's {con['both']}: a bowel-dysfunction definition, not a constipation one.",
        f"Focal and generalized epilepsy share {epi['both']} concepts. Laxative purchases have no condition concept at all.",
    ]
    f.reading(left, y + 6, W - left - 10, lines, strong=(1,))
    f.source("OMOP standard vocabularies from Athena; concept IDs published, SNOMED CT content not redistributed")
    f.save("fig06_snomed_contents")
    paper_single(rows)


def paper_single(rows):
    """The paper's single-column layout: the same squares, one condition per block, no headline or reading."""
    W = sv.SINGLE
    left, right, pitch, sq = 10.0, sv.SINGLE - 8, 11.0, 8.6
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    kx = f.key_row(left, 12, [("square", sv.FIRST, "only first"), ("square", sv.SECOND, "only second")])
    f.polygon([(kx, 6), (kx + 6, 6), (kx, 12)], sv.FIRST)
    f.polygon([(kx + 6, 6), (kx + 6, 12), (kx, 12)], sv.SECOND)
    f.text(kx + 11, 12, "in both", 7.0, sv.INK_2)
    y = 32.0
    for r in rows:
        kinds = ["a"] * r["only_a"] + ["ab"] * r["both"] + ["b"] * r["only_b"]
        assert left + len(kinds) * pitch <= right, "squares run past the column"
        f.text(left, y, r["cond"], 7.6, sv.INK)
        f.text(right, y, f"{r['both']} of {len(kinds)} shared", 7.0, sv.INK_2, anchor="end")
        for k, kind in enumerate(kinds):
            x0, y0 = left + k * pitch, y + 4
            if kind == "ab":
                f.polygon([(x0, y0), (x0 + sq, y0), (x0, y0 + sq)], sv.FIRST, mark=True)
                f.polygon([(x0 + sq, y0), (x0 + sq, y0 + sq), (x0, y0 + sq)], sv.SECOND, mark=True)
            else:
                f.rect(x0, y0, sq, sq, sv.FIRST if kind == "a" else sv.SECOND, rx=1, mark=True)
        n_a, n_b = r["only_a"] + r["both"], r["only_b"] + r["both"]
        f.text(left, y + 22, f"{r['alab']} {n_a}, {r['blab']} {n_b}", 7.0, sv.DIM, max_w=right - left)
        y += 34
        if r["b"] == NOTE_FOR:
            f.text(left, y - 3, "adds irritable bowel syndrome, functional", 7.0, sv.INK_2, max_w=right - left)
            f.text(left, y + 6, "diarrhoea, neurogenic bowel and more", 7.0, sv.INK_2, max_w=right - left)
            y += 18
    f.h = y - 6
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.save("fig06_snomed_contents_paper", paper=False, png=False)


if __name__ == "__main__":
    main()
