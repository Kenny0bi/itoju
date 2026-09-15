"""Figure: sample overlap shows up in the cross-trait LDSC intercept, where it should, and nowhere else.

Each row is one psychiatric GWAS; each dot is its cross-trait intercept with one FinnGen endpoint, laid out
as a beeswarm so no dot hides another. The pale band around zero spans plus or minus 1.96 times the median
intercept standard error: where a GWAS with no people in common with FinnGen should sit. Autism,
schizophrenia, bipolar disorder and clinically defined depression contain no FinnGen samples. The full
depression GWAS and its EHR subset include FinnGen R5, and the PTSD Freeze 3 cohort list includes FinnGen,
so their intercepts should lift off zero. The dark tick is the median.

Source: results/ldsc/rg.tsv (standard two-step LDSC). Every sentence in the reading panel is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv

ROOT = sv.ROOT
BLOCKS = [("NO FINNGEN SAMPLES", ["ASD", "SCZ", "BIP", "MDD_Clin"]),
          ("INCLUDE FINNGEN SAMPLES", ["PTSD", "MDD", "MDD_EHR"])]
NAME = dict(sv.TRAIT_NAME, MDD="depression, all")
LO, HI = -0.02, 0.05
R_DOT = 1.55


def swarm(xs, d):
    placed, out = [], [0.0] * len(xs)
    for i in np.argsort(xs):
        k = 0
        while True:
            off = ((k + 1) // 2) * d * (1 if k % 2 else -1)
            if all((xs[i] - px) ** 2 + (off - py) ** 2 >= d * d for px, py in placed):
                break
            k += 1
        placed.append((xs[i], off))
        out[i] = off
    return out


def main():
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    psych = {"SCZ", "BIP", "MDD", "PTSD", "MDD_EHR", "MDD_Clin", "ASD"}
    endpoints = set(rg[rg.trait1 == "ASD"].trait2) - psych
    sub = rg[rg.trait1.isin(psych) & rg.trait2.isin(endpoints)].dropna(subset=["gcov_int"])
    band = 1.96 * float(sub.gcov_int_se.median())

    W = sv.SINGLE
    left, lab_w = 12.0, 84.0
    ax0, ax1 = left + lab_w, W - 10
    X = sv.scale(LO, HI, ax0, ax1)
    row_h, head_h, top = 26.0, 16.0, 102.0
    body = sum(head_h + row_h * len(ts) for _, ts in BLOCKS)
    H = top + body + 104
    f = sv.Figure(W, H)
    f.header("itoju  /  sample overlap",
             "Shared samples lift the intercept",
             "Cross-trait LDSC intercept with each FinnGen endpoint")
    f.key_row(left, 64, [("band", sv.NOISE, "where no shared samples sit")])
    f.key_row(left, 76, [("dot", sv.INK_2, "one endpoint"), ("line", sv.INK, "median")])
    f.rect(X(-band), top - 4, X(band) - X(-band), body, sv.NOISE)
    f.line(X(0), top - 4, X(0), top + body, sv.RULE, 0.5)

    y = top
    medians = {}
    for head, traits in BLOCKS:
        f.text(left, y + 9, head, 7.0, sv.DIM, spacing=0.4)
        y += head_h
        for t in traits:
            s_ = sub[sub.trait1 == t]
            cy = y + row_h / 2
            xs = np.array([X(np.clip(v, LO, HI)) for v in s_.gcov_int])
            offs = swarm(xs, 2 * R_DOT * 1.04)
            col = sv.TRAIT[t]
            for x, o in zip(xs, offs):
                if t == "MDD_Clin":
                    f.circle(x, cy + o, R_DOT - 0.2, sv.GROUND, stroke=col, sw=0.8)
                else:
                    f.circle(x, cy + o, R_DOT, col)
            med = float(np.median(s_.gcov_int))
            medians[t] = med
            f.line(X(med), cy - 9, X(med), cy + 9, sv.INK, 1.3)
            f.text(left + lab_w - 6, cy + 2.5, NAME[t], 7.0, sv.INK_2, anchor="end")
            y += row_h
    ya = top + body + 2
    f.line(ax0, ya, ax1, ya, sv.RULE, 0.5)
    for v in (-0.02, 0, 0.02, 0.04):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, "0" if v == 0 else sv.fmt(v, 2, sign=True), 7.0, sv.DIM, anchor="middle")
    f.text((ax0 + ax1) / 2, ya + 20, "cross-trait intercept", 7.0, sv.DIM, anchor="middle")

    none = [medians[t] for t in BLOCKS[0][1]]
    print("  medians:", {k: round(v, 4) for k, v in medians.items()}, f"band {band:.4f}")
    lines = [f"Median {min(none):.3f} to {max(none):.3f} with no FinnGen samples;",
             f"{medians['PTSD']:.3f} for PTSD, {medians['MDD']:.3f} and {medians['MDD_EHR']:.3f} for",
             "the depression GWAS that include FinnGen."]
    f.reading(left, ya + 30, W - left - 10, lines, strong=())
    f.source("two-step LDSC, results/ldsc/rg.tsv")
    f.save("fig18_intercepts")


if __name__ == "__main__":
    main()
