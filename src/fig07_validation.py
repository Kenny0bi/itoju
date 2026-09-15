"""Figure: validation. The pipeline lands on published estimates before it is trusted.

Single column, two parts.
a  Each published estimate is drawn the way every comparison in this study is drawn. The pale sleeve,
   centred on the published value, is the range that value allows: its 95% interval, or, when the paper
   prints no SE, the tolerance of plus or minus 0.05 fixed in the analysis plan. The stem runs from the
   published value to this pipeline's estimate, so a stem that stays inside its sleeve reproduces the paper.
   1. Autism SNP heritability, liability scale, K = 0.012 (Grove et al. 2019: 0.118, SE 0.010)
   2. Autism and ADHD genetic correlation (Grove et al. 2019: 0.360; no SE printed)
b  Does Finnish registry autism measure the same condition as the Danish autism GWAS? Each rg is drawn as
   a quantile dot plot of its sampling distribution: 20 dots, each one twentieth of the probability,
   against rg = 1, the value if both measured the same genetics.

Sources: results/ldsc/h2_ASD_Grove2019.log, rg_validation_ASD_ADHD.log, rg.tsv.
"""
import re

import numpy as np
import pandas as pd
from scipy.stats import norm

import itoju_svg as sv

ROOT = sv.ROOT
ASD = sv.TRAIT["ASD"]
PUB_H2, PUB_H2_SE = 0.118, 0.010
PUB_RG, TOL = 0.360, 0.05
BW, NDOT = 0.07, 20
FIN = [("FinnGen autism (F84.0, F84.5)", "KRA_PSY_AUTISM_EXMORE"),
       ("FinnGen autism spectrum (F84)", "KRA_PSY_DEVWIDE_EXMORE")]


def grab(path, label):
    t = (ROOT / "results" / "ldsc" / path).read_text()
    m = re.search(re.escape(label) + r":\s*(-?[\d.]+)\s*\(([\d.]+)\)", t)
    return float(m.group(1)), float(m.group(2))


def main():
    h2, h2se = grab("h2_ASD_Grove2019.log", "Total Liability scale h2")
    rg, rgse = grab("rg_validation_ASD_ADHD.log", "Genetic Correlation")
    assert abs(h2 - PUB_H2) < 1.96 * PUB_H2_SE and abs(rg - PUB_RG) < TOL
    rgt = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    man = pd.read_csv(ROOT / "data" / "raw" / "finngen_R12_manifest.tsv", sep="\t").set_index("phenocode")

    W = sv.SINGLE
    left = 12.0
    ax0, ax1 = 24.0, W - 12
    f = sv.Figure(W, 2000.0)      # drawn tall, trimmed to the content below
    f.header("itoju  /  validation", ["The pipeline reproduces published", "autism estimates first"],
             "checked before any definition was compared")
    f.key_row(left, 76, [("band", sv.NOISE, "range the paper allows"), ("dot", ASD, "this pipeline")])
    f.key_row(left, 88, [("dash", sv.INK_2, "rg = 1: the same genetics")])

    f.text(left, 110, "a  Against the published values", 8.2, sv.INK, family=sv.SERIF)
    X = sv.scale(-0.1, 0.1, ax0, ax1)
    rows = [("Autism SNP heritability, K = 0.012", f"paper {PUB_H2:.3f} (SE {PUB_H2_SE:.3f}); here {h2:.3f}",
             h2 - PUB_H2, 1.96 * PUB_H2_SE),
            ("Autism and ADHD genetic correlation", f"paper {PUB_RG:.3f}, no SE: 0.05 tolerance; here {rg:.3f}",
             rg - PUB_RG, TOL)]
    y = 124.0
    for lab, sub, d, half in rows:
        f.text(left, y, lab, 7.2, sv.INK)
        f.text(left, y + 9.5, sub, 7.0, sv.DIM)
        cy = y + 21
        f.rect(X(-half), cy - 4, X(half) - X(-half), 8, sv.NOISE, rx=1.5)
        f.line(X(0), cy - 6, X(0), cy + 6, sv.RULE, 0.6)
        f.line(X(0), cy, X(d), cy, ASD, 1.6, mark=True)
        f.circle(X(d), cy, 2.8, ASD)
        y += 38
    ya = y - 8
    f.line(ax0, ya, ax1, ya, sv.RULE, 0.5)
    for v, s in ((-0.1, "-0.1"), (-0.05, "-0.05"), (0, "0"), (0.05, "+0.05"), (0.1, "+0.1")):
        f.line(X(v), ya, X(v), ya + 2.5, sv.RULE, 0.5)
        f.text(X(v), ya + 10.5, s, 7.0, sv.DIM, anchor="middle")
    f.text((ax0 + ax1) / 2, ya + 20, "this pipeline minus published", 7.0, sv.DIM, anchor="middle")

    yb = ya + 44
    f.text(left, yb, "b  Finnish registry autism and the GWAS", 8.2, sv.INK, family=sv.SERIF)
    f.text(left, yb + 10, "each dot: 1/20 of the estimate's probability", 7.0, sv.DIM)
    XB = sv.scale(0.2, 1.6, ax0, ax1)
    d = BW / 1.4 * (ax1 - ax0)
    top = yb + 22
    fin_vals = []
    tops = []
    for lab, e in FIN:
        r = rgt[(rgt.trait1 == "ASD") & (rgt.trait2 == e)].iloc[0]
        q = norm.ppf((np.arange(NDOT) + 0.5) / NDOT, r.rg, r.se)
        bins = np.round(q / BW).astype(int)
        stack = int(np.bincount(bins - bins.min()).max())
        base = top + 22 + stack * d
        f.text(left, top + 7, f"{lab}, {int(man.loc[e, 'num_cases']):,} cases", 7.0, sv.INK_2)
        f.text(left, top + 16, f"rg {r.rg:.2f} (SE {r.se:.2f})", 7.0, sv.INK)
        for k in np.unique(bins):
            for j in range(int((bins == k).sum())):
                f.circle(XB(k * BW), base - (j + 0.5) * d, d * 0.43, ASD)
        f.line(ax0, base, ax1, base, sv.GRID, 0.5)
        f.line(XB(1.0), top + 20, XB(1.0), base, sv.INK_2, 0.8, dash="2 1.6", mark=True)
        tops.append(top)
        fin_vals.append((r.rg, r.se))
        top = base + 8
    yx = top - 6
    f.line(ax0, yx, ax1, yx, sv.RULE, 0.5)
    for v in (0.2, 0.6, 1.0, 1.4):
        f.line(XB(v), yx, XB(v), yx + 2.5, sv.RULE, 0.5)
        f.text(XB(v), yx + 10.5, f"{v:g}", 7.0, sv.DIM, anchor="middle")
    f.text((ax0 + ax1) / 2, yx + 20, "rg with iPSYCH-PGC autism", 7.0, sv.DIM, anchor="middle")

    z = [(1 - v) / s for v, s in fin_vals]
    print(f"  h2 {h2:.4f} ({h2se:.4f}) vs {PUB_H2}; rg ADHD {rg:.4f} ({rgse:.4f}) vs {PUB_RG}; "
          f"FinnGen rg {fin_vals}, SE below 1: {[round(v, 2) for v in z]}")
    lines = [f"Heritability {h2:.3f} against {PUB_H2:.3f} published;",
             f"rg with ADHD {rg:.3f} against {PUB_RG:.3f}.",
             f"FinnGen autism sits {z[0]:.1f} SE below rg = 1,",
             f"autism spectrum {z[1]:.1f} SE below it."]
    yr = yx + 32
    end = f.reading(left, yr, W - left - 10, lines, strong=(0,))
    f.h = end + 20
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.source("Grove et al. 2019; results/ldsc")
    f.save("fig07_validation")


if __name__ == "__main__":
    main()
