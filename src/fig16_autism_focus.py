"""Figure: autism, and how large a definition effect these data could have seen.

One row per alternative definition, compared with the first definition of its condition, with the
iPSYCH-PGC autism GWAS as the shared trait. Each row is a caliper on one axis of shift in genetic
correlation: the darker sleeve is plus or minus 1.96 SE of the difference (shifts the test cannot tell
from chance), the paler sleeve is the smallest shift detectable at 80% power, and the violet stem is the
shift itself. The columns at the right give rg under the two definitions (common SNP set, the estimates
the test compares), the shift, and the detectable shift, so a well-powered null reads at a glance.

Source: results/definition_effect/pairwise_delta.tsv. Every sentence in the reading panel is computed.
"""
import pandas as pd

import itoju_svg as sv
from itoju_labels import GROUPS

ROOT = sv.ROOT
SETS = ["Epilepsy", "Sleep apnoea", "Insomnia", "Constipation", "ADHD"]
CAP = 0.6
BONF = 0.05 / 135
ASD = sv.TRAIT["ASD"]


def main():
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw = pw[pw.shared == "ASD"]
    groups = dict(GROUPS)
    rows = []
    for name in SETS:
        defs = groups[name]
        first_e, first_lab = defs[0]
        rows.append(("head", name, first_lab))
        for e, lab in defs[1:]:
            r = pw[((pw.def1 == first_e) & (pw.def2 == e)) | ((pw.def1 == e) & (pw.def2 == first_e))]
            assert len(r) == 1, (first_e, e)
            r = r.iloc[0]
            if r.def1 == first_e:
                rg_a, rg_b, shift = r.rg1, r.rg2, -r.delta
            else:
                rg_a, rg_b, shift = r.rg2, r.rg1, r.delta
            rows.append(("row", name, lab, rg_a, rg_b, shift, 1.96 * r.se_delta, r.min_detectable_delta, r.p))

    W = sv.DOUBLE
    left, ax0, ax1 = 14.0, 182.0, 350.0
    col_rg, col_shift, col_det = 414.0, 450.0, 506.0
    head_h, row_h, top = 16.0, 15.5, 108.0
    body = sum(head_h if r[0] == "head" else row_h for r in rows)
    H = top + body + 108
    f = sv.Figure(W, H)
    f.header("itoju  /  autism",
             "For autism, no change of definition moves the genetics beyond chance",
             "Shift in genetic correlation with iPSYCH-PGC autism against the first definition of each condition")
    f.key_row(left, 72, [("band", sv.NOISE, "chance: 1.96 SE of the difference"),
                         ("band", sv.NOISE_2, "detectable 4 times in 5"),
                         ("dot", ASD, "p < 0.05"), ("hollow", ASD, "not significant")])
    X = sv.scale(-CAP, CAP, ax0, ax1)
    f.text(col_rg, 96, "rg, first to second", 7.0, sv.DIM, anchor="end")
    f.text(col_shift, 96, "shift", 7.0, sv.DIM, anchor="end")
    f.text(col_det, 96, "detectable", 7.0, sv.DIM, anchor="end")

    y = top
    for r in rows:
        if r[0] == "head":
            _, name, first_lab = r
            f.text(left, y + 8, name, 7.6, sv.INK)
            vs = first_lab[0].lower() + first_lab[1:] if not first_lab[:2].isupper() else first_lab
            f.text(left + sv.text_width(name, 7.6) + 8, y + 8, f"vs {vs}", 7.0, sv.DIM)
            y += head_h
            continue
        _, name, lab, rg_a, rg_b, shift, hw, mdd, p = r
        cy = y + 4.5
        for v in (-0.5, 0.5):
            f.line(X(v), y - 3, X(v), y + row_h - 3, sv.GRID, 0.5)
        f.text(left + 10, y + 7, lab, 7.0, sv.INK_2)
        f.rect(X(-min(mdd, CAP)), cy - 3.5, X(min(mdd, CAP)) - X(-min(mdd, CAP)), 7, sv.NOISE_2, rx=1)
        f.rect(X(-min(hw, CAP)), cy - 3.5, X(min(hw, CAP)) - X(-min(hw, CAP)), 7, sv.NOISE, rx=1)
        f.line(X(0), cy - 5.5, X(0), cy + 5.5, sv.RULE, 0.5)
        f.line(X(0), cy, X(max(min(shift, CAP), -CAP)), cy, ASD, 1.6, mark=True)
        if p < 0.05:
            f.circle(X(shift), cy, 2.6, ASD)
        else:
            f.circle(X(shift), cy, 2.3, sv.GROUND, stroke=ASD, sw=1.0)
        f.text(col_rg, y + 7, f"{rg_a:.2f} to {rg_b:.2f}", 7.0, sv.INK_2, anchor="end")
        f.text(col_shift, y + 7, sv.fmt(shift, 2, sign=True), 7.0, sv.INK if p < 0.05 else sv.INK_2, anchor="end")
        f.text(col_det, y + 7, f"{mdd:.2f}", 7.0, sv.INK_2, anchor="end")
        y += row_h
    ya = top + body + 2
    f.line(ax0, ya, ax1, ya, sv.RULE, 0.5)
    for v, s in ((-0.5, "-0.5"), (0, "0"), (0.5, "+0.5")):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, s, 7.0, sv.DIM, anchor="middle")
    f.text((ax0 + ax1) / 2, ya + 20, "shift in rg with autism (sleeves cut at 0.6)", 7.0, sv.DIM, anchor="middle")

    data = [r for r in rows if r[0] == "row"]
    n_nom = sum(r[8] < 0.05 for r in data)
    n_bonf = sum(r[8] < BONF for r in data)
    con = next(r for r in data if r[1] == "Constipation")
    apn = [r for r in data if r[1] == "Sleep apnoea"]
    epi = [r for r in data if r[1] == "Epilepsy"]
    print(f"  autism rows {len(data)}: p<0.05 {n_nom}, Bonferroni {n_bonf}; constipation shift {con[5]:+.4f} detectable {con[7]:.3f}; "
          f"apnoea detectable {min(r[7] for r in apn):.3f}-{max(r[7] for r in apn):.3f}; epilepsy {min(r[7] for r in epi):.3f}-{max(r[7] for r in epi):.3f}")
    lines = [
        f"Against each condition's first definition, {n_nom} of {len(data)} changes reach p < 0.05; {n_bonf} survive Bonferroni.",
        f"Constipation moves by {sv.fmt(con[5], 3, sign=True)} where a shift of {con[7]:.2f} was detectable: a well-powered null.",
        f"Sleep apnoea could show a shift of {min(r[7] for r in apn):.2f}; epilepsy needs {min(r[7] for r in epi):.2f} "
        f"to {max(r[7] for r in epi):.2f}, which these data cannot resolve.",
    ]
    f.reading(left, ya + 30, W - left - 10, lines, strong=(0,))
    f.source("iPSYCH-PGC autism and FinnGen R12; detectable shift = 2.80 SE (alpha 0.05, power 0.80)")
    f.save("fig16_autism_focus")
    paper_single(rows)


def paper_single(rows):
    """The paper's single-column layout: each definition named above its caliper, shift and detectable shift at right."""
    W = sv.SINGLE
    left, right = 10.0, sv.SINGLE - 8
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    f.key_row(left, 12, [("band", sv.NOISE, "chance: 1.96 SE"), ("band", sv.NOISE_2, "detectable 4 in 5")])
    f.key_row(left, 24, [("dot", ASD, "p < 0.05"), ("hollow", ASD, "not significant")])
    X = sv.scale(-CAP, CAP, left + 4, right - 4)
    f.text(right, 42, "shift, detectable", 7.0, sv.DIM, anchor="end")
    top = 56.0
    body_top = top
    y = top
    for r in rows:
        if r[0] == "head":
            _, name, first_lab = r
            vs = first_lab[0].lower() + first_lab[1:] if not first_lab[:2].isupper() else first_lab
            f.text(left, y + 4, name, 7.6, sv.INK)
            f.text(left + sv.text_width(name, 7.6) + 6, y + 4, f"vs {vs}", 7.0, sv.DIM, max_w=right - left - sv.text_width(name, 7.6) - 6)
            y += 12
            continue
        _, name, lab, rg_a, rg_b, shift, hw, mdd, p = r
        f.text(left + 6, y + 4, lab, 7.0, sv.INK_2)
        f.text(right, y + 4, f"{sv.fmt(shift, 2, sign=True)}, {mdd:.2f}", 7.0, sv.INK if p < 0.05 else sv.INK_2, anchor="end")
        cy = y + 11
        f.rect(X(-min(mdd, CAP)), cy - 3, X(min(mdd, CAP)) - X(-min(mdd, CAP)), 6, sv.NOISE_2, rx=1)
        f.rect(X(-min(hw, CAP)), cy - 3, X(min(hw, CAP)) - X(-min(hw, CAP)), 6, sv.NOISE, rx=1)
        f.line(X(0), cy - 4.5, X(0), cy + 4.5, sv.RULE, 0.5)
        f.line(X(0), cy, X(max(min(shift, CAP), -CAP)), cy, ASD, 1.5, mark=True)
        if p < 0.05:
            f.circle(X(shift), cy, 2.4, ASD)
        else:
            f.circle(X(shift), cy, 2.1, sv.GROUND, stroke=ASD, sw=1.0)
        y += 20
    for v in (-0.5, 0.5):
        f.line(X(v), body_top - 2, X(v), y - 2, sv.GRID, 0.4)
    ya = y + 1
    f.line(X(-CAP), ya, X(CAP), ya, sv.RULE, 0.5)
    for v, s in ((-0.5, "-0.5"), (0, "0"), (0.5, "+0.5")):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, s, 7.0, sv.DIM, anchor="middle")
    f.text((left + right) / 2, ya + 20, "shift in rg with autism (sleeves cut at 0.6)", 7.0, sv.DIM, anchor="middle")
    f.h = ya + 26
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.save("fig16_autism_focus_paper", paper=False, png=False)


if __name__ == "__main__":
    main()
