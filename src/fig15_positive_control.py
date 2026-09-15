"""Figure: the positive control. Depression defined from health records casts a wider genetic net.

Depression has two PGC GWAS that differ only in how cases were defined: from electronic health records,
and by clinical assessment. The known result (Cai et al. 2020) is that the looser definition carries less
specific genetics, so it should correlate more strongly with other conditions.

One row per other trait or FinnGen endpoint (31), sorted by the difference. Each arrow starts at the
genetic correlation with clinically defined depression (hollow dot) and ends at the correlation with
EHR-defined depression (arrowhead), so an arrow pointing right means the EHR definition correlates more
strongly. Arrows whose difference reaches p < 0.05 in the shared-block jackknife are drawn heavier in
turmeric. The column at the right sets each difference inside its own noise sleeve (1.96 SE).

Sources: results/definition_effect/pairwise_delta.tsv (family "MDD definition (positive control)"),
results/ldsc/h2.tsv and rg.tsv. Every sentence in the reading panel is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv
from itoju_labels import ENDPOINT_GROUP, ENDPOINT_LABEL, TRAIT_LABEL

ROOT = sv.ROOT
STRONG = "#b98500"       # the depression hue, one step darker so a thin arrow still reads on the ground
WEAK = "#8D9C97"
DCAP = 0.45


def label(e):
    if e in TRAIT_LABEL:
        return TRAIT_LABEL[e]
    g, lab = ENDPOINT_GROUP[e], ENDPOINT_LABEL[e]
    if g in ("Neighbouring conditions", "Constipation") or lab.lower().startswith(g.lower()):
        return lab
    if g == "Intellectual disability":
        return "Intellectual disability (F7)" if "F7)" in lab else "Intellectual disability (F70)"
    return f"{g}: {lab[0].lower() + lab[1:] if not lab[:2].isupper() else lab}"


def arrow(f, x0, x1, y, colour, width):
    f.line(x0, y, x1, y, colour, width, mark=True)
    d = 1 if x1 >= x0 else -1
    if abs(x1 - x0) > 0.5:
        f.polygon([(x1, y), (x1 - d * 4.2, y - 2.4), (x1 - d * 4.2, y + 2.4)], colour, mark=True)


def main():
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    c = pw[pw.family == "MDD definition (positive control)"].copy()
    ehr_first = c.def1 == "MDD_EHR"
    assert set(c.def1) | set(c.def2) == {"MDD_EHR", "MDD_Clin"}
    c["ehr"] = np.where(ehr_first, c.rg1, c.rg2)
    c["clin"] = np.where(ehr_first, c.rg2, c.rg1)
    c["diff"] = c.ehr - c.clin
    c["hw"] = 1.96 * c.se_delta
    c = c.sort_values("diff", ascending=False).reset_index(drop=True)
    n = len(c)
    bonf = 0.05 / n
    h2 = pd.read_csv(ROOT / "results" / "ldsc" / "h2.tsv", sep="\t").set_index("trait")
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    between = rg[((rg.trait1 == "MDD_EHR") & (rg.trait2 == "MDD_Clin")) | ((rg.trait1 == "MDD_Clin") & (rg.trait2 == "MDD_EHR"))].iloc[0]

    W = sv.DOUBLE
    left, lab_r, ax0, ax1 = 14.0, 196.0, 204.0, 386.0
    d0, d1, col_v = 404.0, 478.0, 508.0
    top, row_h = 110.0, 11.0
    H = top + n * row_h + 98
    f = sv.Figure(W, H)
    f.header("itoju  /  the positive control",
             "Looser depression casts a wider genetic net",
             "Genetic correlation with depression defined by clinical assessment and from health records")
    n_up = int((c["diff"] > 0).sum())
    f.text(W - 12, 38, f"{n_up} of {n}", 20.0, sv.INK, family=sv.SERIF, anchor="end")
    f.text(W - 12, 49, "arrows point right", 7.0, sv.DIM, anchor="end")
    kx = f.key_row(left, 74, [("hollow", sv.INK_2, "clinical assessment")])
    f.line(kx, 71.5, kx + 12, 71.5, sv.INK_2, 1.1)
    f.polygon([(kx + 14, 71.5), (kx + 9.8, 69.1), (kx + 9.8, 73.9)], sv.INK_2)
    kx += 18
    f.text(kx, 74, "health records (EHR)", 7.0, sv.INK_2)
    kx += sv.text_width("health records (EHR)") + 10
    f.line(kx, 71.5, kx + 12, 71.5, STRONG, 1.9)
    f.polygon([(kx + 14, 71.5), (kx + 9.8, 69.1), (kx + 9.8, 73.9)], STRONG)
    kx += 18
    f.text(kx, 74, "difference p < 0.05", 7.0, sv.INK_2)
    kx += sv.text_width("difference p < 0.05") + 10
    f.key_row(kx, 74, [("band", sv.NOISE, "noise, 1.96 SE")])

    X = sv.scale(-0.1, 0.9, ax0, ax1)
    D = sv.scale(-DCAP, DCAP, d0, d1)
    f.text((ax0 + ax1) / 2, 98, "genetic correlation with depression", 7.0, sv.DIM, anchor="middle")
    f.text((d0 + d1) / 2, 98, "EHR minus clinical", 7.0, sv.DIM, anchor="middle")
    for v in (0, 0.5):
        f.line(X(v), top - 5, X(v), top + n * row_h, sv.GRID, 0.5)
    for i, r in c.iterrows():
        y = top + i * row_h + 3
        sig = r.p < 0.05
        f.text(lab_r, y + 2.5, label(r.shared), 7.0, sv.INK if sig else sv.INK_2, anchor="end")
        col, wid = (STRONG, 1.9) if sig else (WEAK, 1.1)
        f.circle(X(r.clin), y, 2.2, sv.GROUND, stroke=sv.INK_2 if not sig else STRONG, sw=0.9)
        arrow(f, X(r.clin) + (2.2 if r.ehr >= r.clin else -2.2), X(r.ehr), y, col, wid)
        f.rect(D(-min(r.hw, DCAP)), y - 3, D(min(r.hw, DCAP)) - D(-min(r.hw, DCAP)), 6, sv.NOISE, rx=1)
        f.line(D(0), y - 4.5, D(0), y + 4.5, sv.RULE, 0.5)
        f.circle(D(max(min(r["diff"], DCAP), -DCAP)), y, 2.2, STRONG if sig else sv.GROUND,
                 stroke=None if sig else sv.INK_2, sw=0.9)
        f.text(col_v, y + 2.5, sv.fmt(r["diff"], 2, sign=True), 7.0, sv.INK if sig else sv.DIM, anchor="end")
    ya = top + n * row_h + 2
    f.line(ax0, ya, ax1, ya, sv.RULE, 0.5)
    for v in (0, 0.5):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, f"{v:g}", 7.0, sv.DIM, anchor="middle")
    f.line(d0, ya, d1, ya, sv.RULE, 0.5)
    for v, s in ((-0.4, "-0.4"), (0, "0"), (0.4, "+0.4")):
        f.line(D(v), ya, D(v), ya + 2.5, sv.RULE, 0.5)
        f.text(D(v), ya + 10.5, s, 7.0, sv.DIM, anchor="middle")

    n_nom = int((c.p < 0.05).sum())
    n_bonf = int((c.p < bonf).sum())
    print(f"  control rows {n}: EHR higher {n_up}, p<0.05 {n_nom}, Bonferroni over {n} {n_bonf}; rg between {between.rg:.3f} ({between.se:.3f}); "
          f"h2 clinical {h2.loc['MDD_Clin', 'h2']:.3f}, EHR {h2.loc['MDD_EHR', 'h2']:.3f}; detectable {c.min_detectable_delta.min():.2f}-{c.min_detectable_delta.max():.2f}")
    lines = [
        f"EHR-defined depression correlates more strongly in {n_up} of {n}; {n_nom} reach p < 0.05, {n_bonf} survive Bonferroni.",
        f"The two definitions correlate at {between.rg:.3f} (SE {between.se:.3f}), and clinical heritability is "
        f"{h2.loc['MDD_Clin', 'h2']:.3f} against {h2.loc['MDD_EHR', 'h2']:.3f}.",
        "The two GWAS share almost no people, so their noise does not cancel and single differences stay wide.",
    ]
    f.reading(left, ya + 22, W - left - 10, lines, strong=(0,))
    f.source("PGC MDD2025 EHR and clinical GWAS, FinnGen R12; shared-block jackknife")
    f.save("fig15_positive_control")
    paper_single(c)


def paper_single(c):
    """The paper's single-column layout: the same arrows, shorter labels, no difference column."""
    W = sv.SINGLE
    left, lab_r, ax0, ax1 = 10.0, 126.0, 132.0, W - 8
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    kx = f.key_row(left, 12, [("hollow", sv.INK_2, "clinical assessment")])
    f.line(kx, 9.5, kx + 10, 9.5, sv.INK_2, 1.1)
    f.polygon([(kx + 12, 9.5), (kx + 8, 7.1), (kx + 8, 11.9)], sv.INK_2)
    f.text(kx + 16, 12, "health records (EHR)", 7.0, sv.INK_2)
    f.line(left, 21.5, left + 10, 21.5, STRONG, 1.9)
    f.polygon([(left + 12, 21.5), (left + 8, 19.1), (left + 8, 23.9)], STRONG)
    f.text(left + 16, 24, "difference p < 0.05; right: EHR higher", 7.0, sv.INK_2, max_w=W - 10 - left - 16)
    X = sv.scale(-0.1, 0.9, ax0, ax1)
    top, row_h = 38.0, 10.4
    for v in (0, 0.5):
        f.line(X(v), top - 4, X(v), top + len(c) * row_h, sv.GRID, 0.5)
    for i, r in c.iterrows():
        y = top + i * row_h + 3
        sig = r.p < 0.05
        lab = label(r.shared)
        if sv.text_width(lab) > lab_r - left:          # the group prefix does not fit in one column
            lab = ENDPOINT_LABEL.get(r.shared, lab)
        f.text(lab_r, y + 2.5, lab, 7.0, sv.INK if sig else sv.INK_2, anchor="end", max_w=lab_r - left)
        col, wid = (STRONG, 1.8) if sig else (WEAK, 1.0)
        f.circle(X(r.clin), y, 2.0, sv.GROUND, stroke=STRONG if sig else sv.INK_2, sw=0.9)
        arrow(f, X(r.clin) + (2.0 if r.ehr >= r.clin else -2.0), X(r.ehr), y, col, wid)
    ya = top + len(c) * row_h + 2
    f.line(ax0, ya, ax1, ya, sv.RULE, 0.5)
    for v in (0, 0.5):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, f"{v:g}", 7.0, sv.DIM, anchor="middle")
    f.text((left + ax1) / 2, ya + 20, "genetic correlation with depression", 7.0, sv.DIM, anchor="middle")
    f.h = ya + 26
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.save("fig15_positive_control_paper", paper=False, png=False)


if __name__ == "__main__":
    main()
