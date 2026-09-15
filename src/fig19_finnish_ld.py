"""Figure 19. The LD landscape of chromosomes 21 and 22, seen through two populations.

FinnGen is a Finnish cohort, a population with its own linkage disequilibrium structure, while the
standard LDSC reference is 489 European samples from 1000 Genomes. LD scores built here from the public
FinnGen R12 LD matrix (520,210 Finnish samples; pairs with r2 > 0.01 within 1 Mb, HapMap3 SNPs, plus 1
for the SNP itself) are compared with the standard European scores for the same SNPs.

Single column, stacked.
a, b  Along each chromosome, the mean LD score in 250 kb windows: European above the line, Finnish
      mirrored below it. Where the two shapes mirror each other, the populations tag the genome alike.
c     The per-SNP ratio Finnish / European as a quantile dot plot: 50 dots, each 2% of SNPs, on a log
      axis centered at 1.

Colors validated with the dataviz palette validator (both pass all checks on white).
Sources: data/ref/finngen_ld/{21,22}.l2.ldscore.gz (src/06_finnish_ldscores.py) and
data/ref/z8182036/eur_w_ld_chr/{21,22}.l2.ldscore.gz.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

import viz_style as vs

ROOT = vs.ROOT
FIN = "#2f6fc4"
EUR = "#c9761a"
WIN = 250_000
NDOT = 50
FIG_H = 5.2


def main():
    vs.apply()
    frames = []
    for c in (21, 22):
        fin = pd.read_csv(ROOT / "data" / "ref" / "finngen_ld" / f"{c}.l2.ldscore.gz", sep="\t")
        eur = pd.read_csv(ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr" / f"{c}.l2.ldscore.gz", sep="\t")
        frames.append(fin.merge(eur[["SNP", "L2"]], on="SNP", suffixes=("_fin", "_eur")))
    m = pd.concat(frames)
    m = m[(m.L2_eur > 0) & (m.L2_fin > 0)]
    r = np.corrcoef(m.L2_fin, m.L2_eur)[0, 1]
    ratio = m.L2_fin / m.L2_eur
    med = float(np.median(ratio))
    print(f"  {len(m):,} SNPs, r {r:.3f}, median ratio {med:.3f}, {np.mean(ratio < 1):.1%} lower in Finns")

    fig = plt.figure(figsize=(vs.SINGLE, FIG_H))
    for k, (col, lab) in enumerate(((EUR, "European LD scores (1000 Genomes), above the line"),
                                    (FIN, "Finnish LD scores (FinnGen R12), mirrored below"))):
        y = 0.975 - k * 0.035
        fig.add_artist(Rectangle((0.03, y - 0.011), 0.035, 0.022, transform=fig.transFigure, facecolor=col,
                                 edgecolor="none"))
        fig.text(0.08, y, lab, fontsize=7, color=vs.INK_2, va="center")

    ymax = 0
    wins = {}
    for c in (21, 22):
        s = m[m.CHR == c]
        edges = np.arange(s.BP.min() // WIN * WIN, s.BP.max() + WIN, WIN)
        idx = np.digitize(s.BP, edges)
        centers = (edges[:-1] + WIN / 2) / 1e6
        e_mean = np.full(len(centers), np.nan)
        f_mean = np.full(len(centers), np.nan)
        for k, grp in s.groupby(idx):
            if 1 <= k <= len(centers) and len(grp) >= 3:
                e_mean[k - 1] = grp.L2_eur.mean()
                f_mean[k - 1] = grp.L2_fin.mean()
        wins[c] = (centers, e_mean, f_mean)
        ymax = max(ymax, np.nanmax(e_mean), np.nanmax(f_mean))
    ymax = np.ceil(ymax / 20) * 20

    for i, c in enumerate((21, 22)):
        centers, e_mean, f_mean = wins[c]
        ax = fig.add_axes([0.15, 0.665 - i * 0.28, 0.82, 0.215])
        ok = ~np.isnan(e_mean)
        ax.fill_between(centers, 0, np.where(ok, e_mean, 0), where=ok, color=EUR, lw=0, step="mid", zorder=2)
        ax.fill_between(centers, 0, -np.where(ok, f_mean, 0), where=ok, color=FIN, lw=0, step="mid", zorder=2)
        ax.axhline(0, color="white", lw=0.6, zorder=3)
        ax.set_ylim(-ymax, ymax)
        ticks = [-ymax, 0, ymax]
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{abs(t):g}" for t in ticks])
        ax.set_xlim(np.nanmin(centers) - 0.5, np.nanmax(centers) + 0.5)
        ax.set_ylabel("mean LD score", fontsize=7)
        ax.text(0.01, 0.97, f"{'ab'[i]}  chromosome {c}", transform=ax.transAxes, fontsize=8, color=vs.INK, va="top")
        ax.set_xlabel("position (Mb, GRCh37)", fontsize=7, labelpad=1)

    b = fig.add_axes([0.15, 0.085, 0.44, 0.17])
    lo, hi = -1.5, 1.5
    bw = 0.12
    q = np.log2(np.quantile(ratio, (np.arange(NDOT) + 0.5) / NDOT))
    bins = np.round(q / bw).astype(int)
    xs, hs = [], []
    for k in np.unique(bins):
        n = int((bins == k).sum())
        xs += [k * bw] * n
        hs += list(np.arange(n) + 0.5)
    unit = bw / (hi - lo) * 0.44 * vs.SINGLE * 72
    h_units = 0.17 * FIG_H * 72 / unit
    b.scatter(np.clip(xs, lo, hi), hs, s=(unit * 0.86) ** 2, color=vs.INK_2, linewidths=0, zorder=3)
    b.plot([np.log2(med)] * 2, [0, max(hs) + 1.5], color=FIN, lw=1.2, zorder=2)
    b.set_xlim(lo, hi)
    b.set_ylim(0, h_units)
    b.set_yticks([])
    b.spines["left"].set_visible(False)
    b.set_xticks([-1, 0, 1])
    b.set_xticklabels(["1/2", "1", "2"])
    b.xaxis.itoju_value_of = lambda v: 2.0 ** v
    b.set_xlabel("Finnish / European LD score", fontsize=7)
    fig.text(0.63, 0.25, f"c  SNP by SNP\n{len(m):,} SNPs\nr = {r:.3f}\nblue line: median\nratio {med:.3f}\n"
             "each dot: 2% of SNPs", fontsize=7, color=vs.INK_2, va="top", linespacing=1.3)
    vs.save(fig, "fig19_finnish_ld")


if __name__ == "__main__":
    main()
