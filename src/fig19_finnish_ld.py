"""Figure 19. Finnish LD versus the European reference, on chromosomes 21 and 22.

FinnGen is a Finnish cohort, a population with its own linkage disequilibrium structure, while the
standard LDSC reference is 489 European samples from 1000 Genomes. This compares, SNP by SNP, the
standard European LD scores with LD scores built here from the public FinnGen R12 LD matrix
(520,210 Finnish samples; pairs with r2 > 0.01 within 1 Mb, HapMap3 SNPs, plus 1 for the SNP itself).

a  Log-density hexbin of Finnish against European LD score, with the line of equality.
b  Distribution of the per-SNP ratio Finnish / European, on a log axis centered at 1.

Sources: data/ref/finngen_ld/{21,22}.l2.ldscore.gz (src/06_finnish_ldscores.py) and
data/ref/z8182036/eur_w_ld_chr/{21,22}.l2.ldscore.gz.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import viz_style as vs

ROOT = vs.ROOT


def main():
    vs.apply()
    frames = []
    for c in (21, 22):
        fin = pd.read_csv(ROOT / "data" / "ref" / "finngen_ld" / f"{c}.l2.ldscore.gz", sep="\t")
        eur = pd.read_csv(ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr" / f"{c}.l2.ldscore.gz", sep="\t")
        m = fin.merge(eur[["SNP", "L2"]], on="SNP", suffixes=("_fin", "_eur"))
        m["chr"] = c
        frames.append(m)
    m = pd.concat(frames)
    m = m[(m.L2_eur > 0) & (m.L2_fin > 0)]
    r = np.corrcoef(m.L2_fin, m.L2_eur)[0, 1]
    ratio = m.L2_fin / m.L2_eur

    fig, (a, b) = plt.subplots(1, 2, figsize=(vs.DOUBLE, 2.9), gridspec_kw={"width_ratios": [1, 1.25]})
    fig.subplots_adjust(left=0.08, right=0.98, top=0.86, bottom=0.17, wspace=0.48)
    hb = a.hexbin(m.L2_eur, m.L2_fin, gridsize=55, bins="log", cmap=vs.cmap_seq(), mincnt=1,
                  extent=(0, 120, 0, 120), linewidths=0)
    hb.set_rasterized(True)
    a.plot([0, 120], [0, 120], color=vs.INK_2, lw=0.8)
    a.set_xlim(0, 120)
    a.set_ylim(0, 120)
    a.set_aspect("equal")
    a.set_xlabel("European LD score (1000 Genomes)")
    a.set_ylabel("Finnish LD score (FinnGen R12)")
    n_off = int(((m.L2_eur > 120) | (m.L2_fin > 120)).sum())
    a.text(117, 4, f"{len(m):,} SNPs\nchr 21 and 22\nr = {r:.3f}" + (f"\n{n_off} beyond 120" if n_off else ""),
           fontsize=7, color=vs.INK_2, va="bottom", ha="right")
    a.set_title("a  SNP by SNP (line: equal values)", loc="left", fontsize=8.5)
    cb = fig.colorbar(hb, ax=a, fraction=0.046, pad=0.02)
    cb.set_label("SNPs per hexagon", fontsize=7, color=vs.INK_2)
    cb.ax.tick_params(labelsize=7, colors=vs.MUTED)
    cb.outline.set_visible(False)

    lr = np.log2(ratio.clip(1 / 16, 16))
    bins = np.linspace(-2, 2, 61)
    inside = lr[(lr > -2) & (lr < 2)]
    b.hist(inside, bins=bins, color=vs.SEQ[7], edgecolor="white", linewidth=0.3)
    med = float(np.median(ratio))
    b.axvline(np.log2(med), color=vs.INK, lw=0.9)
    b.set_xticks([-2, -1, 0, 1, 2])
    b.set_xticklabels(["1/4", "1/2", "1", "2", "4"])
    b.xaxis.itoju_value_of = lambda v: 2.0 ** v
    b.set_xlabel("Finnish / European LD score (log scale)")
    b.set_ylabel("SNPs")
    below, above = int((lr <= -2).sum()), int((lr >= 2).sum())
    b.text(0.98, 0.97, f"black line: median ratio {med:.3f}\n{np.mean(ratio < 1):.0%} of SNPs lower in Finns"
           + (f"\n{below} below 1/4 and {above} above 4\nnot shown" if below or above else ""),
           transform=b.transAxes, fontsize=7, color=vs.INK_2, va="top", ha="right")
    b.set_title("b  The ratio", loc="left", fontsize=8.5)
    vs.save(fig, "fig19_finnish_ld")


if __name__ == "__main__":
    main()
