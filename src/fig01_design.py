"""Figure: the study design, walked through once with real numbers.

Four steps, left to right, each drawn with the numbers of one worked example (sleep apnoea and depression),
and each introducing a mark the later figures reuse:
  1  One condition, several borders: the three FinnGen R12 sleep apnoea definitions as bars of cases on one
     scale, ochre for the condition's first definition and blue for the alternatives; the dashed line
     carries the hospital-records length down across the others.
  2  Each border meets five psychiatric GWAS in LD score regression.
  3  The genetic correlation with depression moves when the border moves (common SNP set).
  4  Is the move bigger than its noise? The shift against plus or minus 1.96 SE of the difference from the
     shared-block jackknife: the sleeve of the later figures.
The strip underneath lists the checks that run beside the main design.

Sources: data/raw/finngen_R12_manifest.tsv, results/definition_effect/pairwise_delta.tsv, itoju_labels.
"""
import pandas as pd

import itoju_svg as sv
from itoju_labels import ENDPOINT_ORDER

ROOT = sv.ROOT
APNOEA = [("G6_SLEEPAPNO", "hospital records", "hospital"), ("G6_SLEEPAPNO_INCLAVO", "+ primary care", "+ primary"),
          ("SLEEP", "any sleep disorder", "any sleep")]
MDD = sv.TRAIT["MDD"]
STEP_W, STEP_GAP = 112.0, 13.5


def chevron(f, x, y):
    f.polyline([(x, y - 4), (x + 4, y), (x, y + 4)], sv.RULE, 1.0, mark=False)


def main():
    man = pd.read_csv(ROOT / "data" / "raw" / "finngen_R12_manifest.tsv", sep="\t").set_index("phenocode")
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    row = pw[(pw.shared == "MDD") & (pw.def1 == "G6_SLEEPAPNO") & (pw.def2 == "SLEEP")].iloc[0]
    rg_a, rg_b = float(row.rg1), float(row.rg2)
    shift, hw = rg_b - rg_a, 1.96 * float(row.se_delta)
    assert abs(shift + row.delta) < 1e-9          # delta = rg(def1) - rg(def2)
    n_endpoints = len(ENDPOINT_ORDER)
    print(f"  rg {rg_a:.4f} -> {rg_b:.4f}, shift {shift:+.4f}, 1.96 SE {hw:.4f}, p {row.p:.2e}, endpoints {n_endpoints}")

    W = sv.DOUBLE
    left = 14.0
    f = sv.Figure(W, 350.0)
    f.header("itoju  /  study design",
             "Move a condition's border in the health record, and ask how far the genetics move",
             "One worked example with real numbers: sleep apnoea in FinnGen R12 and depression")
    f.key_row(left, 70, [("square", sv.FIRST, "first definition"), ("square", sv.SECOND, "alternative definitions"),
                         ("band", sv.NOISE, "noise, 1.96 SE of the difference")])
    heads = [("1  Three borders", "sleep apnoea, FinnGen R12", "cases per definition"),
             ("2  Five GWAS each", "LD score regression", f"{n_endpoints} endpoints in all"),
             ("3  The rg moves", "genetic correlation", "with depression"),
             ("4  Beyond its noise?", "shared-block jackknife", "difference test")]
    xs = [left + k * (STEP_W + STEP_GAP) for k in range(4)]
    for x0, (t, a, b) in zip(xs, heads):
        f.text(x0, 98, t, 8.6, sv.INK, family=sv.SERIF)
        f.text(x0, 108.5, a, 7.0, sv.DIM, max_w=STEP_W)
        f.text(x0, 117.5, b, 7.0, sv.DIM, max_w=STEP_W)
    for x0 in xs[:3]:
        chevron(f, x0 + STEP_W + 4, 180)

    # 1 three borders
    x0 = xs[0]
    cases = [int(man.loc[e, "num_cases"]) for e, _, _ in APNOEA]
    k = 78.0 / max(cases)
    for i, ((e, lab, _), n) in enumerate(zip(APNOEA, cases)):
        y = 140 + i * 28
        f.text(x0, y, lab, 7.0, sv.INK)
        f.rect(x0, y + 4, n * k, 8, sv.FIRST if i == 0 else sv.SECOND, rx=1)
        f.text(x0 + n * k + 4, y + 11, f"{n:,}", 7.0, sv.INK_2)
    for i in (1, 2):                                   # only across the bars, never through a label
        yb_ = 140 + i * 28
        f.line(x0 + cases[0] * k, yb_ + 3, x0 + cases[0] * k, yb_ + 15, sv.INK, 0.8, dash="1.6 1.2", mark=True)
    f.text(x0, 224, "dashed: hospital length", 7.0, sv.DIM)

    # 2 each border meets five GWAS
    x0 = xs[1]
    for i, (_, _, short) in enumerate(APNOEA):
        yb = 142 + i * 32
        for j in range(5):
            f.line(x0 + 45, yb + 8, x0 + 54.5, 140 + j * 20, sv.GRID, 0.6)
        f.rect(x0, yb, 45, 16, sv.PANEL_2, rx=2)
        f.rect(x0, yb, 3, 16, sv.FIRST if i == 0 else sv.SECOND)
        f.text(x0 + 6, yb + 10.5, short, 7.0, sv.INK)
    for j, t in enumerate(sv.FIVE):
        f.circle(x0 + 56.5, 140 + j * 20, 2.8, sv.TRAIT[t], mark=False)
        f.text(x0 + 62, 142.5 + j * 20, sv.TRAIT_NAME[t], 7.0, sv.INK_2, max_w=STEP_W - 62 + 6)

    # 3 the correlation moves
    x0 = xs[2]
    X = sv.scale(0.30, 0.46, x0 + 4, x0 + STEP_W - 4)
    ya = 200.0
    f.line(x0, ya, x0 + STEP_W, ya, sv.RULE, 0.6)
    mx = (X(rg_a) + X(rg_b)) / 2
    tip = (X(rg_b), ya - 4.0)
    ctrl = (mx, ya - 58)
    tx, ty = tip[0] - ctrl[0], tip[1] - ctrl[1]
    tn = (tx * tx + ty * ty) ** 0.5
    ux, uy = tx / tn, ty / tn                          # the curve's direction where it lands
    back = (tip[0] - 5.5 * ux, tip[1] - 5.5 * uy)
    f.path(f"M {X(rg_a):.2f},{ya - 4:.2f} Q {ctrl[0]:.2f},{ctrl[1]:.2f} {back[0]:.2f},{back[1]:.2f}", stroke=MDD, width=1.6)
    f.polygon([tip, (back[0] - 2.8 * uy, back[1] + 2.8 * ux), (back[0] + 2.8 * uy, back[1] - 2.8 * ux)], MDD)
    f.circle(X(rg_a), ya, 3.2, sv.FIRST, stroke=sv.GROUND, sw=0.6)
    f.circle(X(rg_b), ya, 3.2, sv.SECOND, stroke=sv.GROUND, sw=0.6)
    f.text(mx, ya - 34, sv.fmt(shift, 3, sign=True), 7.0, sv.INK, anchor="middle")
    f.text(X(rg_a), ya + 12, "hospital", 7.0, sv.INK_2, anchor="middle")
    f.text(X(rg_a), ya + 21, f"rg {rg_a:.3f}", 7.0, sv.INK, anchor="middle")
    f.text(X(rg_b), ya + 12, "any sleep", 7.0, sv.INK_2, anchor="middle")
    f.text(X(rg_b), ya + 21, f"rg {rg_b:.3f}", 7.0, sv.INK, anchor="middle")

    # 4 beyond its noise
    x0 = xs[3]
    D = sv.scale(-0.03, 0.09, x0 + 2, x0 + STEP_W - 2)
    yc = 170.0
    f.text(x0, 146, f"noise: 1.96 SE = {hw:.3f}", 7.0, sv.INK_2, max_w=STEP_W)
    f.rect(D(-hw), yc - 7, D(hw) - D(-hw), 14, sv.NOISE, rx=1.5)
    f.line(D(0), yc - 11, D(0), yc + 11, sv.INK_2, 0.6)
    f.line(D(0), yc, D(shift), yc, MDD, 1.8, mark=True)
    f.circle(D(shift), yc, 3.0, MDD)
    f.text(x0, 198, f"shift {sv.fmt(shift, 3, sign=True)}", 7.0, sv.INK)
    f.text(x0, 207, f"p = {row.p:.1e}", 7.0, sv.INK)
    f.text(x0, 220, "far outside its noise", 7.0, sv.INK_2)

    # checks beside the design
    f.line(left, 242, W - 10, 242, sv.GRID, 0.6)
    f.text(left, 256, "CHECKS BESIDE THE DESIGN", 7.0, sv.DIM, spacing=0.6)
    checks = ["published autism heritability and ADHD rg reproduced", "positive control: EHR against clinical depression",
              "definition codes mapped to SNOMED CT through OMOP", "Finnish against European LD scores"]
    for i, c in enumerate(checks):
        cx = left + (i % 2) * 252
        cy = 268 + (i // 2) * 10
        f.circle(cx + 2.5, cy - 2.5, 2.0, "none", stroke=sv.INK_2, sw=0.8, mark=False)
        f.text(cx + 9, cy, c, 7.0, sv.INK_2, max_w=240)

    lines = [f"Widening sleep apnoea to any sleep disorder moves rg with depression from {rg_a:.3f} to {rg_b:.3f}.",
             f"Pairing the two estimates block by block keeps the noise to 1.96 SE = {hw:.3f}, so p = {row.p:.1e}."]
    f.reading(left, 292, W - left - 10, lines, strong=(0,))
    f.source("FinnGen R12 summary statistics and endpoints; PGC and iPSYCH-PGC GWAS; OMOP vocabularies")
    f.save("fig01_design")
    paper_single(cases, rg_a, rg_b, shift, hw, float(row.p), n_endpoints)


def paper_single(cases, rg_a, rg_b, shift, hw, p, n_endpoints):
    """The paper's single-column layout: the same four steps, two by two, no headline or reading."""
    W = sv.SINGLE
    left = 10.0
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    f.key_row(left, 12, [("square", sv.FIRST, "first definition"), ("square", sv.SECOND, "alternatives")])
    f.key_row(left, 24, [("band", sv.NOISE, "noise, 1.96 SE of the difference")])
    cols = (left, 130.0)
    rows_y = (42.0, 158.0)
    heads = ["1  Three borders", "2  Five GWAS each", "3  The rg moves", "4  Beyond its noise?"]
    for k, h in enumerate(heads):
        f.text(cols[k % 2], rows_y[k // 2], h, 8.0, sv.INK, family=sv.SERIF)

    # 1 the three definitions as bars of cases
    x0, y0 = cols[0], rows_y[0]
    k_ = 58.0 / max(cases)
    for i, ((e, lab, _), n) in enumerate(zip(APNOEA, cases)):
        y = y0 + 14 + i * 17
        f.text(x0, y, lab, 7.0, sv.INK_2)
        f.rect(x0, y + 2.5, n * k_, 6, sv.FIRST if i == 0 else sv.SECOND, rx=1)
        f.text(x0 + n * k_ + 3, y + 8, f"{n:,}", 7.0, sv.INK_2)
    for i in (1, 2):
        yb = y0 + 14 + i * 17
        f.line(x0 + cases[0] * k_, yb + 2, x0 + cases[0] * k_, yb + 9, sv.INK, 0.8, dash="1.6 1.2", mark=True)
    f.text(x0, y0 + 76, "dashed: hospital length", 7.0, sv.DIM)

    # 2 each definition meets five GWAS
    x0, y0 = cols[1], rows_y[0]
    for i, (_, _, short) in enumerate(APNOEA):
        yb = y0 + 14 + i * 20
        for j in range(5):
            f.line(x0 + 38, yb + 7, x0 + 44, y0 + 16 + j * 14, sv.GRID, 0.5)
        f.rect(x0, yb, 38, 14, sv.PANEL_2, rx=2)
        f.rect(x0, yb, 2.5, 14, sv.FIRST if i == 0 else sv.SECOND)
        f.text(x0 + 5, yb + 9.5, short, 7.0, sv.INK)
    for j, t in enumerate(sv.FIVE):
        f.circle(x0 + 46, y0 + 16 + j * 14, 2.5, sv.TRAIT[t], mark=False)
        f.text(x0 + 51, y0 + 18.5 + j * 14, sv.TRAIT_NAME[t], 7.0, sv.INK_2, max_w=W - 10 - (x0 + 51))

    # 3 the correlation moves
    x0, y0 = cols[0], rows_y[1]
    X = sv.scale(0.30, 0.46, x0 + 4, x0 + 102)
    ya = y0 + 58
    f.line(x0, ya, x0 + 108, ya, sv.RULE, 0.6)
    mx = (X(rg_a) + X(rg_b)) / 2
    tip = (X(rg_b), ya - 4.0)
    ctrl = (mx, ya - 46)
    tx, ty = tip[0] - ctrl[0], tip[1] - ctrl[1]
    tn = (tx * tx + ty * ty) ** 0.5
    ux, uy = tx / tn, ty / tn
    back = (tip[0] - 5.0 * ux, tip[1] - 5.0 * uy)
    f.path(f"M {X(rg_a):.2f},{ya - 4:.2f} Q {ctrl[0]:.2f},{ctrl[1]:.2f} {back[0]:.2f},{back[1]:.2f}", stroke=MDD, width=1.5)
    f.polygon([tip, (back[0] - 2.6 * uy, back[1] + 2.6 * ux), (back[0] + 2.6 * uy, back[1] - 2.6 * ux)], MDD)
    f.circle(X(rg_a), ya, 3.0, sv.FIRST, stroke=sv.GROUND, sw=0.6)
    f.circle(X(rg_b), ya, 3.0, sv.SECOND, stroke=sv.GROUND, sw=0.6)
    f.text(mx, ya - 28, sv.fmt(shift, 3, sign=True), 7.0, sv.INK, anchor="middle")
    for x_, name, rg in ((X(rg_a), "hospital", rg_a), (X(rg_b), "any sleep", rg_b)):
        f.text(x_, ya + 11, name, 7.0, sv.INK_2, anchor="middle")
        f.text(x_, ya + 20, f"rg {rg:.3f}", 7.0, sv.INK, anchor="middle")

    # 4 the shift against its noise
    x0, y0 = cols[1], rows_y[1]
    D = sv.scale(-0.03, 0.09, x0 + 2, W - 10)
    yc = y0 + 34
    f.text(x0, y0 + 16, f"noise: 1.96 SE = {hw:.3f}", 7.0, sv.INK_2, max_w=W - 10 - x0)
    f.rect(D(-hw), yc - 6, D(hw) - D(-hw), 12, sv.NOISE, rx=1.5)
    f.line(D(0), yc - 9, D(0), yc + 9, sv.INK_2, 0.6)
    f.line(D(0), yc, D(shift), yc, MDD, 1.6, mark=True)
    f.circle(D(shift), yc, 2.8, MDD)
    f.text(x0, yc + 22, f"shift {sv.fmt(shift, 3, sign=True)}, p = {p:.1e}", 7.0, sv.INK, max_w=W - 10 - x0)
    f.text(x0, yc + 31, "far outside its noise", 7.0, sv.INK_2)

    y = rows_y[1] + 100
    f.line(left, y, W - 10, y, sv.GRID, 0.6)
    f.text(left, y + 12, "CHECKS BESIDE THE DESIGN", 7.0, sv.DIM, spacing=0.6)
    checks = ["published autism heritability and ADHD rg reproduced", "positive control: EHR against clinical depression",
              "definition codes mapped to SNOMED CT through OMOP", "Finnish against European LD scores"]
    for i, c in enumerate(checks):
        cy = y + 24 + i * 9.5
        f.circle(left + 2.5, cy - 2.5, 2.0, "none", stroke=sv.INK_2, sw=0.8, mark=False)
        f.text(left + 9, cy, c, 7.0, sv.INK_2, max_w=W - 10 - left - 9)
    f.h = y + 24 + len(checks) * 9.5
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.save("fig01_design_paper", paper=False, png=False)


if __name__ == "__main__":
    main()
