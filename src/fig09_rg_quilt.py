"""Figure: the genetic correlation quilt, with its headline pattern read first.

a  For each condition with alternative definitions, one line per psychiatric trait: the dot is the median
   genetic correlation across the condition's definitions, the line runs from the lowest to the highest
   definition, so the length of the line is how far the choice of definition alone can move the answer.
b  The full matrix. Each tile is one genetic correlation between a psychiatric GWAS (column) and a FinnGen
   endpoint (row). Colour follows a value-suppressing uncertainty palette: the more precise an estimate,
   the finer its colour steps (8 steps when SE is at most 0.05); as the SE grows, neighbouring values merge
   (4, then 2 steps) and the colour fades toward the noise gray, and above SE 0.20 every value is the same
   gray. The fan key reads outward from uncertain (centre) to precise (outer ring) and left to right from
   negative to positive rg.

Source: results/ldsc/rg.tsv (standard two-step LDSC, full SNP sets). Every sentence in the reading panel is
computed, and the two claims in the headline are asserted.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv
from itoju_labels import GROUPS

ROOT = sv.ROOT
COLS = ["ASD", "SCZ", "BIP", "PTSD", "MDD", "MDD_EHR", "MDD_Clin"]
COL_LABEL = {"ASD": "ASD", "SCZ": "SCZ", "BIP": "BIP", "PTSD": "PTSD", "MDD": "all", "MDD_EHR": "EHR",
             "MDD_Clin": "clin"}
CAP = 0.8
SE_EDGES = [0.05, 0.10, 0.20]          # level 0: SE <= 0.05, 1: <= 0.10, 2: <= 0.20, 3: above
SUPPRESS = [0.0, 0.2, 0.45, 0.88]      # how far each level's colours are pulled toward the noise gray
STOPS = [(0.0, "#1f58b5"), (0.25, "#8db1e6"), (0.5, "#FAF7F1"), (0.625, "#f8c7ab"), (0.75, "#ee8a62"),
         (0.875, "#d44a37"), (1.0, "#971b2c")]
RINGS = [(44.0, 55.0), (33.6, 44.0), (23.2, 33.6), (12.8, 23.2)]
A0, A1 = 158.0, 22.0
ROW, HEAD, GROUP_GAP, PITCH = 9.5, 11.0, 4.0, 17.0


def rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], float)


def hexc(a):
    return "#" + "".join(f"{int(round(v)):02x}" for v in np.clip(a, 0, 255))


def ramp(t):
    for (t0, c0), (t1, c1) in zip(STOPS, STOPS[1:]):
        if t <= t1:
            u = (t - t0) / (t1 - t0)
            return rgb(c0) * (1 - u) + rgb(c1) * u
    return rgb(STOPS[-1][1])


def level(se):
    return int(np.searchsorted(SE_EDGES, se, side="left"))


def vsup(v, lev):
    nb = 8 >> lev
    i = min(nb - 1, int((np.clip(v, -CAP, CAP) + CAP) / (2 * CAP) * nb))
    vc = -CAP + (i + 0.5) * 2 * CAP / nb
    return hexc(ramp((vc + CAP) / (2 * CAP)) * (1 - SUPPRESS[lev]) + rgb(sv.NOISE) * SUPPRESS[lev])


def wedge(cx, cy, r0, r1, a_lo, a_hi):
    p = lambda r, a: (cx + r * np.cos(np.radians(a)), cy - r * np.sin(np.radians(a)))
    (x1, y1), (x2, y2), (x3, y3), (x4, y4) = p(r1, a_hi), p(r1, a_lo), p(r0, a_lo), p(r0, a_hi)
    return (f"M {x1:.2f},{y1:.2f} A {r1} {r1} 0 0 1 {x2:.2f},{y2:.2f} L {x3:.2f},{y3:.2f} "
            f"A {r0} {r0} 0 0 0 {x4:.2f},{y4:.2f} Z")


def main():
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    rg = rg[rg.trait1.isin(COLS)]
    cell = {(r.trait1, r.trait2): r for r in rg.itertuples()}
    conds = [(g, [e for e, _ in items]) for g, items in GROUPS if g != "Neighbouring conditions"]
    five = rg[rg.trait1.isin(sv.FIVE)].dropna(subset=["rg"])
    stat = {}
    for g, es in conds:
        for t in sv.FIVE:
            v = five[(five.trait1 == t) & five.trait2.isin(es)].rg
            stat[(g, t)] = (float(v.median()), float(v.min()), float(v.max()))
    top = {t: max(conds, key=lambda c: stat[(c[0], t)][0])[0] for t in sv.FIVE}
    low = {t: min(conds, key=lambda c: stat[(c[0], t)][0])[0] for t in sv.FIVE}
    n_top = sum(top[t] in ("ADHD", "Insomnia") for t in sv.FIVE)
    n_low = sum(low[t] in ("Epilepsy", "Sleep apnoea") for t in sv.FIVE)
    asd_epi = five[(five.trait1 == "ASD") & five.trait2.isin(dict(conds)["Epilepsy"])]
    widest = max(stat, key=lambda k: stat[k][2] - stat[k][1])
    print("  strongest:", top, "\n  weakest:", low, f"\n  autism x epilepsy max rg {asd_epi.rg.max():.3f} over {len(asd_epi)}",
          f"\n  widest spread: {widest} {stat[widest][1]:.3f} to {stat[widest][2]:.3f}")
    assert n_top >= 4 and n_low == 5 and asd_epi.rg.max() < 0

    W = sv.DOUBLE
    left = 14.0
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    f.header("itoju  /  the whole matrix",
             "ADHD and insomnia share the most genetics; epilepsy and sleep apnoea the least",
             "Genetic correlation of every FinnGen endpoint with every psychiatric GWAS, full SNP sets")
    f.key_traits(left, 70)

    # ---- a: each condition, across its definitions
    ay = 98.0
    f.text(left, ay, "a  Each condition, across its definitions", 8.2, sv.INK, family=sv.SERIF)
    f.text(left, ay + 10, "dot: median rg; line: lowest to highest definition", 7.0, sv.DIM)
    px0, px1 = 104.0, 244.0
    X = sv.scale(-0.4, 0.9, px0, px1)
    sub_h, gap = 6.5, 8.0
    y0 = ay + 24
    y_end = y0 + len(conds) * (5 * sub_h + gap) - gap
    for v in (0.4, 0.8):
        f.line(X(v), y0 - 3, X(v), y_end + 2, sv.GRID, 0.5)
    f.line(X(0), y0 - 3, X(0), y_end + 2, sv.RULE, 0.6)
    y = y0
    for i, (g, _) in enumerate(conds):
        words = [g] if sv.text_width(g) <= px0 - left - 8 else g.split(" ", 1)
        mid = y + 2.5 * sub_h
        for k, w_ in enumerate(words):
            f.text(left, mid + 2.5 + (k - (len(words) - 1) / 2) * 9, w_, 7.0, sv.INK_2)
        for k, t in enumerate(sv.FIVE):
            med, lo, hi = stat[(g, t)]
            cy = y + (k + 0.5) * sub_h
            f.line(X(lo), cy, X(hi), cy, sv.TRAIT[t], 1.2, cap="round", mark=True, opacity=0.6)
            f.circle(X(med), cy, 2.3, sv.TRAIT[t], stroke=sv.GROUND, sw=0.5)
        y += 5 * sub_h + gap
        if i < len(conds) - 1:
            f.line(left, y - gap / 2, px1, y - gap / 2, sv.GRID, 0.4)
    ya = y_end + 4
    f.line(px0, ya, px1, ya, sv.RULE, 0.5)
    for v in (-0.4, 0, 0.4, 0.8):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, "0" if v == 0 else f"{v:+.1f}", 7.0, sv.DIM, anchor="middle")
    f.text((px0 + px1) / 2, ya + 20, "rg with the psychiatric trait", 7.0, sv.DIM, anchor="middle")

    # ---- the colour key for b, under a
    ky = ya + 44
    f.text(left, ky, "TILE COLOUR IN b", 7.0, sv.DIM, spacing=0.6)
    cx, cy = 92.0, ky + 76
    for lev, (r0, r1) in enumerate(RINGS):
        nb = 8 >> lev
        for i in range(nb):
            a_hi = A0 - i * (A0 - A1) / nb
            a_lo = A0 - (i + 1) * (A0 - A1) / nb
            vc = -CAP + (i + 0.5) * 2 * CAP / nb
            f.path(wedge(cx, cy, r0, r1, a_lo, a_hi), stroke=sv.GROUND, width=0.8, fill=vsup(vc, lev))
    for v, lab, anchor in ((-0.8, "-0.8", "end"), (0.0, "0", "middle"), (0.8, "0.8+", "start")):
        a = np.radians(A0 - (v + CAP) / (2 * CAP) * (A0 - A1))
        off = 3 if anchor == "middle" else 0
        f.text(cx + 59 * np.cos(a), cy - 59 * np.sin(a) + 2.4 - off, lab, 7.0, sv.INK_2, anchor=anchor)
    klines = ["hue: rg, blue below zero, red above",
              "colour steps: 8 to SE 0.05, 4 to 0.10, 2 to 0.20",
              "above SE 0.20 every value is one gray"]
    for i, s in enumerate(klines):
        f.text(left, cy + 18 + i * 9, s, 7.0, sv.INK_2, max_w=px1 - left)
    a_bottom = cy + 18 + 2 * 9

    # ---- b: every endpoint
    bx = 262.0
    lab_r, tx0 = 380.0, 384.0
    f.text(bx, ay, "b  Every endpoint", 8.2, sv.INK, family=sv.SERIF)
    f.text(bx, ay + 10, "one tile per genetic correlation", 7.0, sv.DIM)
    d0, d1 = tx0 + COLS.index("MDD") * PITCH + 1, tx0 + (COLS.index("MDD_Clin") + 1) * PITCH - 1
    f.polyline([(d0, ay + 26), (d0, ay + 23), (d1, ay + 23), (d1, ay + 26)], sv.RULE, 0.6, mark=False)
    f.text((d0 + d1) / 2, ay + 20, "depression", 7.0, sv.INK_2, anchor="middle")
    for j, t in enumerate(COLS):
        xc = tx0 + (j + 0.5) * PITCH
        f.text(xc, ay + 35, COL_LABEL[t], 7.0, sv.INK_2, anchor="middle")
        f.circle(xc, ay + 42, 2.4, sv.TRAIT[t], mark=False)
    y = ay + 50
    counts = [0, 0, 0, 0]
    for gname, items in GROUPS:
        f.text(lab_r, y + 8, gname, 7.8, sv.INK, family=sv.SERIF, anchor="end")
        y += HEAD
        for e, lab in items:
            missing = 0
            for j, t in enumerate(COLS):
                r = cell.get((t, e))
                if r is None or pd.isna(r.rg):
                    missing += 1
                    continue
                lev = level(r.se)
                counts[lev] += 1
                f.rect(tx0 + j * PITCH + 0.8, y + 0.8, PITCH - 1.6, ROW - 1.6, vsup(r.rg, lev), rx=0.8, mark=True)
            f.text(lab_r, y + 7, lab, 7.0, sv.INK_2, anchor="end", max_w=lab_r - bx)
            if missing == len(COLS):
                f.text(tx0 + len(COLS) * PITCH / 2, y + 7, "no estimate: h2 below 0", 7.0, sv.DIM, anchor="middle")
            y += ROW
        y += GROUP_GAP
    b_bottom = y

    print("  tiles per SE level:", counts)
    lines = [
        f"Strongest: ADHD or insomnia for {n_top} of 5 psychiatric traits. Weakest: epilepsy or sleep apnoea for all {n_low}.",
        f"Autism and epilepsy correlate below zero under all {len(asd_epi)} epilepsy definitions (highest {asd_epi.rg.max():.2f}).",
        f"{counts[0]} of {sum(counts)} tiles are precise enough (SE to 0.05) for all 8 colour steps; {counts[3]} fade to gray.",
    ]
    end = f.reading(left, max(a_bottom, b_bottom) + 10, W - left - 10, lines, strong=(0,))
    f.h = end + 20
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.source("results/ldsc/rg.tsv, two-step LDSC on full SNP sets; value-suppressing uncertainty palette")
    f.save("fig09_rg_quilt")


if __name__ == "__main__":
    main()
