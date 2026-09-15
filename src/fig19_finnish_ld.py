"""Figure: the LD landscape of chromosomes 21 and 22, seen through two populations.

FinnGen is a Finnish cohort, a population with its own linkage disequilibrium structure, while the standard
LDSC reference is 489 European samples from 1000 Genomes. LD scores built here from the public FinnGen R12
LD matrix (520,210 Finnish samples; pairs with r2 > 0.01 within 1 Mb, HapMap3 SNPs, plus 1 for the SNP
itself) are compared with the standard European scores for the same SNPs.

Single column, stacked.
a, b  Along each chromosome, the mean LD score in 250 kb windows: European above the line, Finnish mirrored
      below it. Where the two shapes mirror each other, the populations tag the genome alike.
c     The per-SNP ratio Finnish / European as a quantile dot plot: 50 dots, each 2% of SNPs, on a log axis
      centred at 1.

Sources: data/ref/finngen_ld/{21,22}.l2.ldscore.gz (src/06_finnish_ldscores.py) and
data/ref/z8182036/eur_w_ld_chr/{21,22}.l2.ldscore.gz. Every sentence in the reading panel is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv

ROOT = sv.ROOT
FIN = sv.SECOND
EUR = sv.FIRST
WIN = 250_000
NDOT = 50
HALF = 34.0


def main():
    frames = []
    for c in (21, 22):
        fin = pd.read_csv(ROOT / "data" / "ref" / "finngen_ld" / f"{c}.l2.ldscore.gz", sep="\t")
        eur = pd.read_csv(ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr" / f"{c}.l2.ldscore.gz", sep="\t")
        frames.append(fin.merge(eur[["SNP", "L2"]], on="SNP", suffixes=("_fin", "_eur")))
    m = pd.concat(frames)
    m = m[(m.L2_eur > 0) & (m.L2_fin > 0)]
    r = float(np.corrcoef(m.L2_fin, m.L2_eur)[0, 1])
    ratio = m.L2_fin / m.L2_eur
    med = float(np.median(ratio))
    lower = float(np.mean(ratio < 1))
    print(f"  {len(m):,} SNPs, r {r:.3f}, median ratio {med:.3f}, {lower:.1%} lower in Finns")

    wins, ymax = {}, 0.0
    for c in (21, 22):
        s = m[m.CHR == c]
        edges = np.arange(s.BP.min() // WIN * WIN, s.BP.max() + WIN, WIN)
        idx = np.digitize(s.BP, edges)
        e_mean = np.zeros(len(edges) - 1)
        f_mean = np.zeros(len(edges) - 1)
        for k, grp in s.groupby(idx):
            if 1 <= k <= len(e_mean) and len(grp) >= 3:
                e_mean[k - 1] = grp.L2_eur.mean()
                f_mean[k - 1] = grp.L2_fin.mean()
        wins[c] = (edges / 1e6, e_mean, f_mean)
        ymax = max(ymax, e_mean.max(), f_mean.max())
    ymax = float(np.ceil(ymax / 20) * 20)

    W = sv.SINGLE
    left = 12.0
    ax0, ax1 = 34.0, W - 10
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    title = (["Finnish and European LD scores", "rise and fall together"] if r > 0.9
             else ["Finnish and European LD scores", "part ways"])
    f.header("itoju  /  Finnish LD", title, "same SNPs, two reference populations")
    f.key_row(left, 76, [("band", EUR, "European, 1000 Genomes: above")])
    f.key_row(left, 88, [("band", FIN, "Finnish, FinnGen R12: mirrored below")])

    t0 = 112.0
    for i, c in enumerate((21, 22)):
        edges, e_mean, f_mean = wins[c]
        f.text(left, t0, f"{'ab'[i]}  Chromosome {c}", 8.2, sv.INK, family=sv.SERIF)
        f.text(left + sv.text_width(f"{'ab'[i]}  Chromosome {c}", 8.2, serif=True) + 6, t0,
               "mean LD score, 250 kb windows", 7.0, sv.DIM)
        y0 = t0 + 12 + HALF
        X = sv.scale(edges[0] - 0.5, edges[-1] + 0.5, ax0, ax1)
        k = HALF / ymax
        for sign, vals, col in ((-1, e_mean, EUR), (1, f_mean, FIN)):
            pts = [(X(edges[0]), y0)]
            for a, b, v in zip(edges[:-1], edges[1:], vals):
                pts += [(X(a), y0 + sign * v * k), (X(b), y0 + sign * v * k)]
            pts.append((X(edges[-1]), y0))
            f.polygon(pts, col, mark=True)
        f.line(ax0, y0, ax1, y0, sv.GROUND, 0.7)
        for v, yy in ((ymax, y0 - HALF), (0, y0), (ymax, y0 + HALF)):
            f.line(ax0 - 2.5, yy, ax0, yy, sv.RULE, 0.5)
            f.text(ax0 - 4, yy + 2.4, f"{v:g}", 7.0, sv.DIM, anchor="end")
        yb = y0 + HALF + 7
        f.line(ax0, yb, ax1, yb, sv.RULE, 0.5)
        for v in range(10, 60, 10):
            if edges[0] - 0.5 <= v <= edges[-1] + 0.5:
                f.line(X(v), yb, X(v), yb + 2.5, sv.RULE, 0.5)
                f.text(X(v), yb + 10.5, f"{v}", 7.0, sv.DIM, anchor="middle")
        f.text(ax1, yb + 20, "position, Mb (GRCh37)", 7.0, sv.DIM, anchor="end")
        t0 = yb + 40

    f.text(left, t0, "c  SNP by SNP", 8.2, sv.INK, family=sv.SERIF)
    f.text(left + sv.text_width("c  SNP by SNP", 8.2, serif=True) + 6, t0, "each dot: 2% of SNPs", 7.0, sv.DIM)
    lo, hi, bw = -1.5, 1.5, 0.12
    XB = sv.scale(lo, hi, ax0, ax1)
    d = bw / (hi - lo) * (ax1 - ax0)
    q = np.log2(np.quantile(ratio, (np.arange(NDOT) + 0.5) / NDOT))
    assert q.min() > lo and q.max() < hi
    bins = np.round(q / bw).astype(int)
    stack = int(np.bincount(bins - bins.min()).max())
    base = t0 + 12 + stack * d + 4
    f.line(XB(np.log2(med)), base - stack * d - 4, XB(np.log2(med)), base, FIN, 1.2, mark=True)
    f.text(XB(np.log2(med)) + 4, base - stack * d - 1, f"median {med:.3f}", 7.0, FIN)
    for kk in np.unique(bins):
        for j in range(int((bins == kk).sum())):
            f.circle(XB(kk * bw), base - (j + 0.5) * d, d * 0.43, sv.INK_2)
    f.line(ax0, base, ax1, base, sv.RULE, 0.5)
    for v, s in ((-1, "1/2"), (0, "1"), (1, "2")):
        f.line(XB(v), base, XB(v), base + 2.5, sv.RULE, 0.5)
        f.text(XB(v), base + 10.5, s, 7.0, sv.DIM, anchor="middle")
    f.text((ax0 + ax1) / 2, base + 20, "Finnish / European LD score, log scale", 7.0, sv.DIM, anchor="middle")

    lines = [f"{len(m):,} SNPs; r = {r:.3f} between panels.",
             f"Median Finnish / European ratio {med:.3f};",
             f"{lower:.0%} of SNPs score lower in Finns."]
    end = f.reading(left, base + 32, W - left - 10, lines, strong=(0,))
    f.h = end + 20
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.source("FinnGen R12 LD matrix; 1000 Genomes EUR")
    f.save("fig19_finnish_ld")


if __name__ == "__main__":
    main()
