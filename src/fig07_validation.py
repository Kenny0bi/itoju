"""Figure 7. Validation: the pipeline reproduces published estimates before it is trusted.

Each row is one estimate: the value this pipeline produced (filled, with a 95% interval) against
the value printed in the source paper (hollow, with its 95% interval where the paper gives an SE).
  1. Autism SNP heritability, liability scale, K = 0.012 (Grove et al. 2019: 0.118, SE 0.010)
  2. Autism x ADHD genetic correlation (Grove et al. 2019: 0.360; no SE printed in the text)
A second panel checks that the Finnish registry agrees with the Danish autism GWAS: rg between
iPSYCH-PGC autism and the two FinnGen autism endpoints, which should be high if both measure the
same condition.

Sources: results/ldsc/h2_ASD_Grove2019.log, rg_validation_ASD_ADHD.log, rg.tsv.
"""
import re

import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs

ROOT = vs.ROOT
Z = 1.96


def grab(path, label):
    t = (ROOT / "results" / "ldsc" / path).read_text()
    m = re.search(re.escape(label) + r":\s*(-?[\d.]+)\s*\(([\d.]+)\)", t)
    return float(m.group(1)), float(m.group(2))


def main():
    vs.apply()
    h2, h2se = grab("h2_ASD_Grove2019.log", "Total Liability scale h2")
    rg, rgse = grab("rg_validation_ASD_ADHD.log", "Genetic Correlation")
    rows = [("Autism SNP heritability\n(liability, K = 0.012)", h2, h2se, 0.118, 0.010),
            ("Autism x ADHD\ngenetic correlation", rg, rgse, 0.360, None)]
    rgt = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    fin = [("FinnGen autism\n(F84.0, F84.5; 888 cases)", "KRA_PSY_AUTISM_EXMORE"),
           ("FinnGen autism spectrum\n(F84; 1,179 cases)", "KRA_PSY_DEVWIDE_EXMORE")]

    fig, (a, b) = plt.subplots(1, 2, figsize=(vs.DOUBLE, 1.9), gridspec_kw={"width_ratios": [1, 1]})
    fig.subplots_adjust(left=0.19, right=0.98, top=0.80, bottom=0.25, wspace=0.75)
    for i, (lab, est, se, pub, pubse) in enumerate(rows):
        y = len(rows) - 1 - i
        a.plot([est - Z * se, est + Z * se], [y + 0.12, y + 0.12], color=vs.TRAIT["ASD"], lw=1.4)
        a.scatter(est, y + 0.12, s=30, color=vs.TRAIT["ASD"], zorder=3)
        if pubse:
            a.plot([pub - Z * pubse, pub + Z * pubse], [y - 0.12, y - 0.12], color=vs.MUTED, lw=1.4)
        a.scatter(pub, y - 0.12, s=30, facecolor="white", edgecolor=vs.INK_2, linewidths=1.0, zorder=3)
        a.text(-0.02, y, lab, transform=a.get_yaxis_transform(), ha="right", va="center", fontsize=7.5, color=vs.INK_2)
        a.text(max(est + Z * se, pub) + 0.012, y, f"{est:.3f} vs {pub:.3f}", va="center", fontsize=7, color=vs.INK_2)
    a.set_yticks([])
    a.spines["left"].set_visible(False)
    a.set_xlim(0, 0.62)
    a.set_ylim(-0.6, len(rows) - 0.4)
    a.set_xlabel("estimate with 95% interval")
    a.set_title("a  Reproducing the autism GWAS paper", loc="left", fontsize=8.5, pad=14)
    a.scatter(0.0, 1.10, s=22, color=vs.TRAIT["ASD"], transform=a.transAxes, clip_on=False)
    a.text(0.03, 1.10, "this pipeline", transform=a.transAxes, fontsize=7, color=vs.INK_2, va="center")
    a.scatter(0.40, 1.10, s=22, facecolor="white", edgecolor=vs.INK_2, transform=a.transAxes, clip_on=False)
    a.text(0.43, 1.10, "published", transform=a.transAxes, fontsize=7, color=vs.INK_2, va="center")

    for i, (lab, e) in enumerate(fin):
        r = rgt[(rgt.trait1 == "ASD") & (rgt.trait2 == e)].iloc[0]
        y = len(fin) - 1 - i
        b.plot([r.rg - Z * r.se, min(r.rg + Z * r.se, 1.4)], [y, y], color=vs.TRAIT["ASD"], lw=1.4)
        b.scatter(r.rg, y, s=30, color=vs.TRAIT["ASD"], zorder=3)
        b.text(-0.02, y, lab, transform=b.get_yaxis_transform(), ha="right", va="center", fontsize=7.5, color=vs.INK_2)
        b.text(0.02, y + 0.30, f"rg {r.rg:.2f} (SE {r.se:.2f})", ha="left", fontsize=7, color=vs.INK_2)
    b.axvline(1, color=vs.MUTED, lw=0.6)
    b.set_yticks([])
    b.spines["left"].set_visible(False)
    b.set_xlim(0, 1.45)
    b.set_ylim(-0.6, len(fin) - 0.2)
    b.set_xlabel("rg with iPSYCH-PGC autism (95% interval)")
    b.set_title("b  Finnish registry autism agrees", loc="left", fontsize=8.5, pad=14)
    vs.save(fig, "fig07_validation")


if __name__ == "__main__":
    main()
