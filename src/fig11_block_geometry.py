"""Figure 11. Why a difference between two genetic correlations can be precise.

Top row: the 200 genome blocks of LDSC's jackknife, each a dot at its rg pseudovalue under the two
definitions (centered on the estimates). Bottom row: a quantile dot plot of the per-block
differences, 20 dots, each standing for 10 blocks. Above the line, the blocks paired as they really
are; below it, the same blocks with the pairing broken (every block of one definition against every
block of the other), which is what the difference would look like if the two estimates shared no
people. The narrower the upper stack, the smaller the SE of the difference.

Pairs: depression x sleep apnoea (nested definitions), PTSD x constipation (partly shared),
PTSD x EHR vs clinical depression (almost nothing shared). Source: results/definition_effect/ldsc.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import viz_style as vs
from itoju_labels import ENDPOINT_LABEL


ROOT = vs.ROOT
RUNS = ROOT / "results" / "definition_effect" / "ldsc"


def pseudo(first, second, run_first, tag=""):
    log = RUNS / f"rg_{run_first}{tag}.log"
    text = log.read_text()
    lines = text[text.find("Summary of Genetic Correlation Results"):].splitlines()[1:]
    hdr = lines[0].split()
    for l in lines[1:]:
        rec = dict(zip(hdr, l.split()))
        if rec["p2"].endswith(f"/{second}.sumstats.gz"):
            break
    stem = f"{log.stem}{run_first}.sumstats.gz_{second}.sumstats.gz"
    g = np.loadtxt(RUNS / f"{stem}.gencov.delete").reshape(-1)
    h1 = np.loadtxt(RUNS / f"{stem}.hsq1.delete").reshape(-1)
    h2 = np.loadtxt(RUNS / f"{stem}.hsq2.delete").reshape(-1)
    rg, se = float(rec["rg"]), float(rec["se"])
    n = len(g)
    return n * rg - (n - 1) * g / np.sqrt(h1 * h2), rg, se


PANELS = [
    ("a", "Nested definitions", "depression x sleep apnoea", "MDD", "G6_SLEEPAPNO", "SLEEP", "MDD", ""),
    ("b", "Partly shared", "PTSD x constipation", "PTSD", "K11_CONSTIPATION", "K11_OTHFUNC", "PTSD", ""),
    ("c", "Almost nothing shared", "PTSD x depression", "PTSD", "MDD_EHR", "MDD_Clin", "PTSD", "__control"),
]
LABEL = dict(ENDPOINT_LABEL, MDD_EHR="EHR depression", MDD_Clin="clinical depression",
             G6_SLEEPAPNO="sleep apnoea, hospital", SLEEP="any sleep disorder")
LIM = 3.2
NDOT = 20
BW = 0.25
FIG_H = 4.9
COLW, GAP, LEFT = 0.262, 0.058, 0.075


def stacks(values):
    q = np.quantile(values, (np.arange(NDOT) + 0.5) / NDOT)
    bins = np.round(q / BW).astype(int)
    xs, hs = [], []
    for b in np.unique(bins):
        k = int((bins == b).sum())
        xs += [b * BW] * k
        hs += list(np.arange(k) + 0.5)
    return np.array(xs), np.array(hs)


def main():
    vs.apply()
    rows = []
    for letter, head, sub, shared, d1, d2, run, tag in PANELS:
        p1, r1, _ = pseudo(shared, d1, run, tag)
        p2, r2, _ = pseudo(shared, d2, run, tag)
        x, y = p1 - r1, p2 - r2
        n = len(x)
        rows.append(dict(
            letter=letter, head=head, sub=sub, color=vs.TRAIT[shared], d1=d1, d2=d2, x=x, y=y,
            corr=np.corrcoef(x, y)[0, 1],
            se_pair=np.sqrt(np.var(p2 - p1, ddof=1) / n),
            se_unp=np.sqrt((np.var(p1, ddof=1) + np.var(p2, ddof=1)) / n),
            paired=stacks(y - x), unpaired=stacks((y[None, :] - x[:, None]).ravel())))
    top_units = max(r["paired"][1].max() for r in rows) + 1.0
    bot_units = max(r["unpaired"][1].max() for r in rows) + 1.0

    fig = plt.figure(figsize=(vs.DOUBLE, FIG_H))
    w_pt = COLW * vs.DOUBLE * 72
    unit_pt = BW / (2 * LIM) * w_pt
    bh = (top_units + bot_units) * unit_pt / (FIG_H * 72)     # dot units map to points 1:1 in x and y
    by = 0.075
    ty = by + bh + 0.165
    th = COLW * vs.DOUBLE / FIG_H
    for i, r in enumerate(rows):
        x0 = LEFT + i * (COLW + GAP)
        top = fig.add_axes([x0, ty, COLW, th])
        t = np.array([-LIM, LIM])
        top.plot(t, t, color=vs.INK_2, lw=0.6, zorder=2)
        inside = (np.abs(r["x"]) <= LIM) & (np.abs(r["y"]) <= LIM)
        top.scatter(r["x"][inside], r["y"][inside], s=7, color=r["color"], alpha=0.85, linewidths=0, zorder=3)
        top.set_xlim(-LIM, LIM)
        top.set_ylim(-LIM, LIM)
        top.set_aspect("equal")
        top.set_xticks([-3, 0, 3])
        top.set_yticks([-3, 0, 3])
        top.set_xlabel(LABEL[r["d1"]], fontsize=7.5, labelpad=2)
        top.set_ylabel(LABEL[r["d2"]], fontsize=7.5, labelpad=1)
        top.text(0.0, 1.15, f"{r['letter']}  {r['head']}", transform=top.transAxes, fontsize=8.5, color=vs.INK)
        top.text(0.0, 1.05, r["sub"], transform=top.transAxes, fontsize=7.5, color=vs.INK_2)
        top.text(0.04, 0.96, f"blocks agree, r = {r['corr']:.2f}", transform=top.transAxes, fontsize=7,
                 color=vs.INK_2, va="top")
        off = int((~inside).sum())
        if off:
            top.text(0.96, 0.04, f"{off} block off scale", transform=top.transAxes, fontsize=7, color=vs.MUTED,
                     ha="right", va="bottom")

        bot = fig.add_axes([x0, by, COLW, bh])
        d_pt = unit_pt * 0.86
        px, ph = r["paired"]
        ux, uh = r["unpaired"]
        bot.scatter(px, ph, s=d_pt ** 2, color=r["color"], linewidths=0, zorder=3)
        bot.scatter(ux, -uh, s=(d_pt * 0.9) ** 2, facecolor="white", edgecolor=vs.MUTED, linewidths=0.7, zorder=3)
        bot.axhline(0, color=vs.INK_2, lw=0.6)
        bot.set_xlim(-LIM, LIM)
        bot.set_ylim(-bot_units, top_units)
        bot.set_yticks([])
        bot.spines["left"].set_visible(False)
        bot.set_xticks([-3, 0, 3])
        bot.set_xlabel("per-block difference", fontsize=7.5, labelpad=2)
        bot.text(-LIM, top_units - 0.3, f"paired\nSE {r['se_pair']:.3f}", fontsize=7, color=vs.INK_2, va="top",
                 linespacing=1.2)
        bot.text(-LIM, -bot_units + 0.3, f"pairing broken\nSE {r['se_unp']:.3f}", fontsize=7, color=vs.INK_2,
                 va="bottom", linespacing=1.2)
    fig.text(LEFT, by + bh + 0.03, "Each dot below stands for 10 of the 200 genome blocks.  Filled: the blocks "
             "paired as they are.\nHollow: every block of one definition against every block of the other, "
             "as if the two estimates shared no people.\nDot colors follow the psychiatric trait: depression in a, PTSD in b and c.",
             fontsize=7, color=vs.INK_2, linespacing=1.3)
    vs.save(fig, "fig11_block_geometry")


if __name__ == "__main__":
    main()
