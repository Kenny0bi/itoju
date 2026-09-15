"""Figure: what moves the genetics when a definition changes, and what does not.

Rows are six changes to a FinnGen definition, grouped by the kind of change (the data source, the
scope of diagnoses let in, the strictness of the rule). Columns are the five psychiatric GWAS. In each
cell the gray sleeve spans plus or minus 1.96 SE of the difference from the shared-block jackknife,
centred on zero: a shift inside it cannot be told from chance. The stem is the shift itself (rg under
the second definition minus rg under the first), in the trait's colour, with its value underneath.
A filled dot has p < 0.05; a ringed dot also survives Bonferroni over all 135 comparisons.

Source: results/definition_effect/pairwise_delta.tsv. Every sentence in the reading panel is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv

ROOT = sv.ROOT
BONF = 0.05 / 135
ROWS = [  # condition, change, kind, first definition, second definition
    ("Sleep apnoea", "add primary care records", "source", "G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO"),
    ("Sleep apnoea", "widen to any sleep disorder", "scope", "G6_SLEEPAPNO", "SLEEP"),
    ("Insomnia", "widen F51.0, G47.0 to all of F51", "scope", "F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"),
    ("Constipation", "widen K59.0 or laxatives to K59", "scope", "K11_CONSTIPATION", "K11_OTHFUNC"),
    ("ADHD", "widen F90.0 to all of F90", "scope", "F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"),
    ("Focal epilepsy", "exclude anyone also generalized", "strictness", "FE", "FE_STRICT"),
]
KIND_LABEL = {"source": "CHANGES THE SOURCE", "scope": "CHANGES THE SCOPE", "strictness": "CHANGES THE RULE"}
LO, HI = -0.12, 0.26


def shift_row(pw, trait, a, b):
    r = pw[(pw.shared == trait) & (((pw.def1 == a) & (pw.def2 == b)) | ((pw.def1 == b) & (pw.def2 == a)))]
    assert len(r) == 1, (trait, a, b, len(r))
    r = r.iloc[0]
    shift = -r.delta if r.def1 == a else r.delta          # delta = rg(def1) - rg(def2)
    return float(shift), float(1.96 * r.se_delta), float(r.p)


def main():
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    cells = {(i, t): shift_row(pw, t, a, b) for i, (_, _, _, a, b) in enumerate(ROWS) for t in sv.FIVE}

    W = sv.DOUBLE
    label_w, col_w, left = 150.0, 67.0, 14.0
    top_rows, row_h = 112.0, 36.0
    H = top_rows + len(ROWS) * row_h + 98
    f = sv.Figure(W, H)
    f.header("itoju  /  what moves the genetics",
             "Widening a definition moves the genetics. Changing its source does not.",
             "Shift in genetic correlation when only the FinnGen definition changes, five psychiatric GWAS")
    f.key_row(left, 66, [("band", sv.NOISE, "noise: 1.96 SE of the difference"),
                         ("dot", sv.INK_2, "p < 0.05"),
                         ("hollow", sv.INK_2, "not significant")])
    # ringed dot drawn by hand so its key matches the mark exactly
    kx = sv.text_width("noise: 1.96 SE of the difference") + sv.text_width("p < 0.05") + sv.text_width("not significant") + 3 * 21 + left + 1
    f.circle(kx + 3, 63.5, 2.8, sv.INK_2, mark=False)
    f.circle(kx + 3, 63.5, 4.6, "none", stroke=sv.INK, sw=0.7, mark=False)
    f.text(kx + 11, 66, "survives Bonferroni (135 tests)", 7.0, sv.INK_2)

    x_cols = [left + label_w + k * col_w for k in range(len(sv.FIVE))]
    for x, t in zip(x_cols, sv.FIVE):
        cx = x + (col_w - 10) / 2 + 2
        wname = sv.text_width(sv.TRAIT_NAME[t])
        f.circle(cx - (wname + 9) / 2 + 2.8, 90, 2.8, sv.TRAIT[t], mark=False)
        f.text(cx - (wname + 9) / 2 + 9, 92.5, sv.TRAIT_NAME[t], 7.0, sv.INK)
        f.keyed = True
    for i, (cond, change, kind, a, b) in enumerate(ROWS):
        y = top_rows + i * row_h
        if i == 0 or kind != ROWS[i - 1][2]:
            f.line(left, y - 4, W - 10, y - 4, sv.RULE, 0.5)
        f.text(left, y + 8, cond, 7.6, sv.INK)
        f.text(left, y + 17.5, change, 7.0, sv.INK_2)
        f.text(left, y + 27, KIND_LABEL[kind], 7.0, sv.DIM, spacing=0.6)
        for x, t in zip(x_cols, sv.FIVE):
            s = sv.scale(LO, HI, x + 4, x + col_w - 6)
            shift, hw, p = cells[(i, t)]
            cy = y + 11
            f.rect(s(max(-hw, LO)), cy - 5, s(min(hw, HI)) - s(max(-hw, LO)), 10, sv.NOISE, rx=1.5)
            f.line(s(0), cy - 8, s(0), cy + 8, sv.RULE, 0.5)
            col = sv.TRAIT[t]
            f.line(s(0), cy, s(shift), cy, col, 2.0, mark=True)
            if p < 0.05:
                f.circle(s(shift), cy, 2.8, col)
            else:
                f.circle(s(shift), cy, 2.4, sv.GROUND, stroke=col, sw=1.0)
            if p < BONF:
                f.circle(s(shift), cy, 4.6, "none", stroke=sv.INK, sw=0.7)
            f.text(s(shift), y + 27, sv.fmt(shift, 3, sign=True), 7.0, sv.INK if p < 0.05 else sv.DIM, anchor="middle")
    ybot = top_rows + len(ROWS) * row_h - 2
    for x in x_cols:
        s = sv.scale(LO, HI, x + 4, x + col_w - 6)
        f.line(x + 4, ybot, x + col_w - 6, ybot, sv.RULE, 0.5)
        for v, lab in ((0, "0"), (0.2, "+0.2")):
            f.line(s(v), ybot, s(v), ybot + 2.5, sv.RULE, 0.5)
            f.text(s(v), ybot + 10, lab, 7.0, sv.DIM, anchor="middle")

    # the reading, computed
    scope = [cells[(i, t)] for i, r in enumerate(ROWS) if r[2] == "scope" for t in sv.FIVE]
    source = [cells[(0, t)] for t in sv.FIVE]
    strict = [cells[(5, t)] for t in sv.FIVE]
    n_scope_pos = sum(c[0] > 0 for c in scope)
    n_scope_nom = sum(c[2] < 0.05 for c in scope)
    n_scope_bonf = sum(c[2] < BONF for c in scope)
    src_max = max(abs(c[0]) for c in source)
    src_hw_min = min(c[1] for c in source)
    src_nom = sum(c[2] < 0.05 for c in source)
    src_bonf = sum(c[2] < BONF for c in source)
    strict_nom = sum(c[2] < 0.05 for c in strict)
    asd_bonf = sum(cells[(i, "ASD")][2] < BONF for i in range(len(ROWS)))
    print(f"  scope cells {len(scope)}: positive {n_scope_pos}, p<0.05 {n_scope_nom}, Bonferroni {n_scope_bonf}; "
          f"source max |shift| {src_max:.4f}, narrowest noise {src_hw_min:.4f}; strict p<0.05 {strict_nom}; autism Bonferroni {asd_bonf}")
    lines = [
        f"Widening the scope: {n_scope_pos} of {len(scope)} shifts are positive, {n_scope_nom} reach p < 0.05, "
        f"{n_scope_bonf} survive Bonferroni.",
        f"Adding primary care: every shift is {src_max:.3f} or smaller; {src_nom} of 5 reach p < 0.05, "
        f"{src_bonf} survive Bonferroni.",
        f"Tightening the focal epilepsy rule: {strict_nom} of 5 shifts reach p < 0.05.",
        f"Autism: {asd_bonf} of 6 survive correction.",
    ]
    yr = top_rows + len(ROWS) * row_h + 16
    f.reading(left, yr, W - left - 10, lines, strong=(0,))
    f.source("FinnGen R12 and PGC summary statistics; shared-block jackknife on 953,474 SNPs")
    f.save("fig10_what_moves")


if __name__ == "__main__":
    main()
