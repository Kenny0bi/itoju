"""Figure: how much do two definitions overlap? It depends on the vocabulary you look through.

For every pair of alternative definitions of the same condition, overlap is measured four ways, from the
codes a registry stores to the clinical meaning behind them:
  ICD-10 codes         the WHO ICD-10 codes each definition includes
  phecodeX groups      the research groupings those codes fall into
  SNOMED CT concepts   the standard concepts the codes map to in the OMOP vocabularies
  SNOMED CT + parents  the same concepts widened to their ancestors within two levels
Each overlap is drawn as an eclipse: two equal discs (ochre the first definition, blue the second) placed so
that the dark shared area, divided by the area the two discs cover together, equals the Jaccard index. Apart
means nothing shared; one dark disc means identical. The geometry is checked against the value it stands
for. Pairs with identical profiles across all four lenses share a row, and the last column counts them.

Source: results/omop/pairwise_omop.tsv (src/08b). Every sentence in the reading panel is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv

ROOT = sv.ROOT
SETS = {
    "Epilepsy": ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"],
    "Sleep apnoea": ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"],
    "Insomnia": ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"],
    "Constipation": ["K11_CONSTIPATION", "K11_OTHFUNC"],
    "ADHD": ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"],
    "Intellectual disability": ["F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE"],
    "Bladder": ["N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD"],
}
SHORT = {"G6_EPLEPSY": "any", "FE": "focal", "FE_STRICT": "focal strict", "FE_MODE": "focal mode",
         "GE": "generalized", "GE_STRICT": "generalized strict", "GE_MODE": "generalized mode",
         "G6_SLEEPAPNO": "hospital", "G6_SLEEPAPNO_INCLAVO": "+ primary care", "SLEEP": "any sleep disorder",
         "F5_INSOMNIA": "F51.0 or G47.0", "KRA_PSY_SLEEP_NONORG_EXMORE": "all of F51",
         "K11_CONSTIPATION": "K59.0 or laxatives", "K11_OTHFUNC": "all of K59", "F5_ADHD": "F90.0",
         "KRA_PSY_HYPERKIN_EXMORE": "all of F90", "F5_MILDRET": "mild (F70)", "KRA_PSY_MENTALRET_EXMORE": "any (F7)",
         "N14_NEUROMUSCDYSBLADD": "neurogenic (N31)", "N14_OTHBLADD": "other (N32)"}
LENSES = [("icd10_jaccard", "ICD-10", "codes"), ("phecode_jaccard", "phecodeX", "groups"),
          ("snomed_jaccard", "SNOMED CT", "concepts"), ("snomed_hier_jaccard", "SNOMED CT", "+ parents")]
R = 4.4
SHARED = sv.INK
PITCH = 12.5
C0, COLW = 236.0, 62.0


def lens_area(d, r=R):
    if d >= 2 * r:
        return 0.0
    return 2 * r * r * np.arccos(d / (2 * r)) - (d / 2) * np.sqrt(4 * r * r - d * d)


def separation(j, r=R):
    """Distance between centres of two equal discs whose shared area over covered area equals j."""
    if j <= 0:
        return 2.25 * r
    if j >= 1:
        return 0.0
    lo, hi = 0.0, 2 * r
    for _ in range(60):
        mid = (lo + hi) / 2
        a = lens_area(mid, r)
        if a / (2 * np.pi * r * r - a) > j:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def eclipse(f, cx, cy, j):
    d = separation(j)
    f.circle(cx - d / 2, cy, R, sv.FIRST, stroke=sv.GROUND, sw=0.5)
    f.circle(cx + d / 2, cy, R, sv.SECOND, stroke=sv.GROUND, sw=0.5)
    if j >= 1:
        f.circle(cx, cy, R, SHARED, stroke=sv.GROUND, sw=0.5)
    elif d < 2 * R:
        h = np.sqrt(R * R - (d / 2) ** 2)
        f.path(f"M {cx:.3f},{cy - h:.3f} A {R} {R} 0 0 1 {cx:.3f},{cy + h:.3f} A {R} {R} 0 0 1 {cx:.3f},{cy - h:.3f} Z",
               stroke="none", width=0, fill=SHARED)
        a = lens_area(d)
        assert abs(a / (2 * np.pi * R * R - a) - j) < 1e-3     # the drawn geometry equals the value


def main():
    pw = pd.read_csv(ROOT / "results" / "omop" / "pairwise_omop.tsv", sep="\t")
    keys = [k for k, _, _ in LENSES]
    rows, allpairs = [], []
    for fam, defs in SETS.items():
        sub = pw[pw.def1.isin(defs) & pw.def2.isin(defs)].copy()
        allpairs.append(sub)
        sub[keys] = sub[keys].round(3)
        for vals, g in sub.groupby(keys, sort=False):
            first = g.iloc[0]
            rows.append((fam, f"{SHORT[first.def1]} vs {SHORT[first.def2]}", list(vals), len(g)))
    ap = pd.concat(allpairs)
    n_pairs = len(ap)
    means = [float(ap[k].mean()) for k in keys]
    up = int((ap.snomed_hier_jaccard > ap.icd10_jaccard + 1e-9).sum())
    down = int((ap.snomed_hier_jaccard < ap.icd10_jaccard - 1e-9).sum())
    same = int(((ap.icd10_jaccard == 0) & (ap.snomed_hier_jaccard > 0)).sum())
    print(f"  pairs {n_pairs}, rows {len(rows)}, means {[round(m, 3) for m in means]}, hier>icd {up}, hier<icd {down}, "
          f"no ICD-10 code shared but SNOMED+parents overlap {same}")

    W = sv.DOUBLE
    left = 14.0
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    title = ("Definitions look more alike through SNOMED CT than through code lists" if means[3] > means[0]
             else "Definitions look less alike through SNOMED CT than through code lists")
    f.header("itoju  /  three lenses", title,
             "Overlap between two definitions of one condition (Jaccard), through four vocabularies")
    f.key_row(left, 70, [("dot", sv.FIRST, "first definition"), ("dot", sv.SECOND, "second definition"),
                         ("dot", SHARED, "what they share, over everything either includes")])
    x = left
    for v, lab in ((0.0, "0: nothing shared"), (1 / 3, "0.33: a third"), (1.0, "1: identical")):
        eclipse(f, x + 10, 83, v)
        f.text(x + 23, 85.5, lab, 7.0, sv.INK_2)
        x += 23 + sv.text_width(lab) + 16

    hy = 110.0
    for k, (_, a, b) in enumerate(LENSES):
        f.text(C0 + k * COLW, hy, a, 7.0, sv.INK)
        f.text(C0 + k * COLW, hy + 9, b, 7.0, sv.INK)
    f.text(W - 10, hy + 9, "pairs", 7.0, sv.DIM, anchor="end")
    f.line(left, hy + 14, W - 10, hy + 14, sv.RULE, 0.5)

    y, last = hy + 26, None
    for fam, label, vals, n in rows:
        if fam != last:
            if last is not None:
                y += 3
            f.text(left, y + 2.5, fam, 7.8, sv.INK, family=sv.SERIF)
            y += PITCH
            last = fam
        f.text(left + 10, y + 2.5, label, 7.0, sv.INK_2, max_w=C0 - left - 20)
        for k, v in enumerate(vals):
            cx = C0 + k * COLW + 10
            eclipse(f, cx, y, v)
            f.text(cx + 14, y + 2.5, f"{v:.2f}", 7.0, sv.INK_2)
        f.text(W - 10, y + 2.5, str(n), 7.0, sv.INK_2 if n > 1 else sv.DIM, anchor="end")
        y += PITCH

    lines = [
        f"Mean overlap over {n_pairs} pairs: ICD-10 {means[0]:.2f}, phecodeX {means[1]:.2f}, SNOMED CT {means[2]:.2f}, "
        f"with parents {means[3]:.2f}.",
        f"Adding SNOMED CT parents raises the overlap above ICD-10 in {up} of {n_pairs} pairs and lowers it in {down}.",
        f"{same} pairs share no ICD-10 code at all, yet overlap once SNOMED CT parents are added.",
    ]
    end = f.reading(left, y + 4, W - left - 10, lines, strong=(0,))
    f.h = end + 20
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.source("results/omop/pairwise_omop.tsv; OMOP vocabularies from Athena")
    f.save("fig05_three_lenses")


if __name__ == "__main__":
    main()
