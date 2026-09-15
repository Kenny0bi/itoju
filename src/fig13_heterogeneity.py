"""Figure: where the definition matters, with each p-value drawn as the tail it is.

Rows are the five testable definition sets, columns the five psychiatric GWAS. Each cell asks whether the
genetic correlation with that trait varies across the definitions of that condition more than the
jackknife covariance allows (Cochran's Q with the full covariance matrix, df = definitions - 1). The curve
is the chi-square distribution Q would follow if the definitions made no difference; the tick is the
observed Q, and the tinted tail beyond it is the p-value, drawn as the area it is. A dot on the tick marks
p < 0.05; an arrow at the right edge marks a Q beyond the drawn range.

Source: results/definition_effect/heterogeneity.tsv. Every sentence in the reading panel is computed.
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2

import itoju_svg as sv

ROOT = sv.ROOT
SETS = [("Sleep apnoea", 3), ("Insomnia", 2), ("Constipation", 2), ("ADHD", 2), ("Epilepsy", 7)]


def main():
    h = pd.read_csv(ROOT / "results" / "definition_effect" / "heterogeneity.tsv", sep="\t")
    cell = {(r.shared, r.family): r for r in h.itertuples()}

    W = sv.DOUBLE
    left, lab_w = 14.0, 104.0
    cw = (W - left - 10 - lab_w) / 5
    top, rh, ch = 112.0, 62.0, 34.0
    H = top + len(SETS) * rh + 80
    f = sv.Figure(W, H)
    f.header("itoju  /  where the definition matters",
             "Definitions disagree for sleep, insomnia and epilepsy; never for ADHD or autism",
             "Cochran's Q across all definitions of one condition, with the full jackknife covariance")
    f.key_row(left, 68, [("line", sv.INK_2, "where Q falls if the definitions made no difference"),
                         ("band", "#E6C9DA", "the p-value: the tail beyond the observed Q")])
    kx = f.key_row(left, 81, [("line", sv.INK, "observed Q, p of 0.05 or more"),
                              ("dot", sv.TRAIT["PTSD"], "observed Q, p < 0.05, in the trait's colour")])
    f.line(kx, 78.5, kx + 10, 78.5, sv.TRAIT["PTSD"], 1.3)
    f.polygon([(kx + 14, 78.5), (kx + 10, 76.1), (kx + 10, 80.9)], sv.TRAIT["PTSD"])
    f.text(kx + 18, 81, "Q beyond the drawn range", 7.0, sv.INK_2)
    for j, t in enumerate(sv.FIVE):
        cx = left + lab_w + j * cw + cw / 2
        wname = sv.text_width(sv.TRAIT_NAME[t])
        f.circle(cx - (wname + 9) / 2 + 2.8, 95.5, 2.8, sv.TRAIT[t], mark=False)
        f.text(cx - (wname + 9) / 2 + 9, 98, sv.TRAIT_NAME[t], 7.0, sv.INK)

    sig = {}
    for i, (fam, k) in enumerate(SETS):
        df = k - 1
        xmax = chi2.ppf(0.995, df) * 1.15
        xs = np.linspace(xmax / 400, xmax, 400)
        dens = chi2.pdf(xs, df)
        ymax = chi2.pdf(0.35, df) * 1.1 if df == 1 else dens.max() * 1.15
        y = top + i * rh
        base = y + ch + 6
        f.text(left, y + 20, fam, 7.6, sv.INK)
        f.text(left, y + 30, f"{k} definitions, df {df}", 7.0, sv.DIM)
        if i:
            f.line(left, y - 4, W - 10, y - 4, sv.GRID, 0.5)
        for j, t in enumerate(sv.FIVE):
            r = cell[(t, fam)]
            x0 = left + lab_w + j * cw + 6
            wpx = cw - 16
            X = sv.scale(0, xmax, x0, x0 + wpx)
            Y = lambda d, base=base, ymax=ymax: base - min(d, ymax) / ymax * ch
            col = sv.TRAIT[t]
            clear = r.p < 0.05
            sig.setdefault(fam, []).append((t, clear, r.Q, r.p))
            if r.Q < xmax:
                m = xs >= r.Q
                tail = [(X(r.Q), base)] + [(X(a), Y(b)) for a, b in zip(xs[m], dens[m])] + [(X(xmax), base)]
                f.polygon(tail, col, opacity=0.32)
            inside = dens <= ymax
            f.polyline([(X(a), Y(b)) for a, b in zip(xs[inside], dens[inside])], sv.INK_2, 0.8, mark=False)
            f.line(x0, base, x0 + wpx, base, sv.RULE, 0.5)
            if r.Q < xmax:
                f.line(X(r.Q), base, X(r.Q), base - ch * 0.85, col if clear else sv.INK, 1.3 if clear else 0.8)
                if clear:
                    f.circle(X(r.Q), base - ch * 0.85, 2.3, col)
            else:
                ya = base - ch * 0.45
                f.line(x0 + wpx * 0.55, ya, x0 + wpx - 3, ya, col, 1.3)
                f.polygon([(x0 + wpx, ya), (x0 + wpx - 4, ya - 2.4), (x0 + wpx - 4, ya + 2.4)], col)
            ptxt = "p<0.001" if r.p < 0.001 else f"p {r.p:.3f}"
            f.text(x0, base + 10, f"Q {r.Q:.1f}", 7.0, sv.INK if clear else sv.DIM)
            f.text(x0 + wpx + 4, base + 10, ptxt, 7.0, sv.INK if clear else sv.DIM, anchor="end")

    n = {fam: sum(c for _, c, _, _ in v) for fam, v in sig.items()}
    epi_traits = [sv.TRAIT_NAME[t] for t, c, _, _ in sig["Epilepsy"] if c]
    asd = sum(c for fam in sig for t, c, _, _ in sig[fam] if t == "ASD")
    print("  p<0.05 per set:", n, "| epilepsy:", epi_traits, "| autism:", asd)
    lines = [f"Sleep apnoea: definitions disagree for {n['Sleep apnoea']} of 5 traits; insomnia {n['Insomnia']} of 5; "
             f"constipation {n['Constipation']} of 5.",
             f"Epilepsy: {n['Epilepsy']} of 5 ({', '.join(epi_traits)}). ADHD: {n['ADHD']} of 5.",
             f"Autism: definitions disagree for {asd} of the 5 conditions."]
    f.reading(left, top + len(SETS) * rh + 6, W - left - 10, lines, strong=(0,))
    f.source("results/definition_effect/heterogeneity.tsv")
    f.save("fig13_heterogeneity")
    paper_single(cell)


def draw_cell(f, r, t, df, x0, wpx, base, ch):
    """One chi-square cell: the null curve, the tail beyond the observed Q, and the tick or arrow."""
    xmax = chi2.ppf(0.995, df) * 1.15
    xs = np.linspace(xmax / 400, xmax, 400)
    dens = chi2.pdf(xs, df)
    ymax = chi2.pdf(0.35, df) * 1.1 if df == 1 else dens.max() * 1.15
    X = sv.scale(0, xmax, x0, x0 + wpx)
    Y = lambda d: base - min(d, ymax) / ymax * ch
    col = sv.TRAIT[t]
    clear = r.p < 0.05
    if r.Q < xmax:
        m = xs >= r.Q
        tail = [(X(r.Q), base)] + [(X(a), Y(b)) for a, b in zip(xs[m], dens[m])] + [(X(xmax), base)]
        f.polygon(tail, col, opacity=0.32)
    inside = dens <= ymax
    f.polyline([(X(a), Y(b)) for a, b in zip(xs[inside], dens[inside])], sv.INK_2, 0.8, mark=False)
    f.line(x0, base, x0 + wpx, base, sv.RULE, 0.5)
    if r.Q < xmax:
        f.line(X(r.Q), base, X(r.Q), base - ch * 0.85, col if clear else sv.INK, 1.3 if clear else 0.8)
        if clear:
            f.circle(X(r.Q), base - ch * 0.85, 2.3, col)
    else:
        ya = base - ch * 0.45
        f.line(x0 + wpx * 0.45, ya, x0 + wpx - 3, ya, col, 1.3)
        f.polygon([(x0 + wpx, ya), (x0 + wpx - 4, ya - 2.4), (x0 + wpx - 4, ya + 2.4)], col)
    return clear


def paper_single(cell):
    """The paper's single-column layout: the same cells, the set named above each row, no headline or reading."""
    W = sv.SINGLE
    left, right = 10.0, sv.SINGLE - 8
    cw = (right - left) / 5
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    f.key_row(left, 12, [("line", sv.INK_2, "Q if definitions made no difference")])
    f.key_row(left, 24, [("band", "#E6C9DA", "p-value: the tail beyond Q")])
    kx = f.key_row(left, 36, [("line", sv.INK, "Q, p of 0.05 or more"), ("dot", sv.TRAIT["PTSD"], "p < 0.05")])
    f.line(left - 1, 45.5, left + 8, 45.5, sv.TRAIT["PTSD"], 1.3)
    f.polygon([(left + 11, 45.5), (left + 7, 43.1), (left + 7, 47.9)], sv.TRAIT["PTSD"])
    f.text(left + 15, 48, "Q beyond the drawn range; colour: trait", 7.0, sv.INK_2)
    hy = 66.0
    for j, t in enumerate(sv.FIVE):
        cx = left + j * cw + cw / 2
        lab = {"ASD": "ASD", "SCZ": "SCZ", "BIP": "BIP", "MDD": "MDD", "PTSD": "PTSD"}[t]
        wl = sv.text_width(lab)
        f.circle(cx - (wl + 8) / 2 + 2.6, hy - 2.5, 2.6, sv.TRAIT[t], mark=False)
        f.text(cx - (wl + 8) / 2 + 8, hy, lab, 7.0, sv.INK)
    y = hy + 16
    for i, (fam, k) in enumerate(SETS):
        f.text(left, y, f"{fam}, {k} definitions", 7.4, sv.INK, family=sv.SERIF)
        base = y + 32
        for j, t in enumerate(sv.FIVE):
            r = cell[(t, fam)]
            x0 = left + j * cw + 3
            clear = draw_cell(f, r, t, k - 1, x0, cw - 8, base, 24)
            ink = sv.INK if clear else sv.DIM
            f.text(x0, base + 9, f"Q {r.Q:.1f}", 7.0, ink, max_w=cw - 2)
            f.text(x0, base + 18, "p<0.001" if r.p < 0.001 else f"p {r.p:.3f}", 7.0, ink, max_w=cw - 2)
        y = base + 32
        if i < len(SETS) - 1:
            f.line(left, y - 9, right, y - 9, sv.GRID, 0.5)
    f.h = y - 8
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.save("fig13_heterogeneity_paper", paper=False, png=False)


if __name__ == "__main__":
    main()
