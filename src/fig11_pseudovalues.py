"""Figure 11. Why a difference between two genetic correlations can be precise.

LDSC's 200-block jackknife produces one pseudovalue of rg per block. On a common SNP set, block k is
the same stretch of genome for both definitions, so each dot here is one block: its rg under
definition 1 (x) against its rg under definition 2 (y). When the two definitions share most of
their cases, the dots lie along the diagonal and most sampling noise cancels in the difference.
The ribbon shows the spread of the per-block difference; its width sets SE(delta).

Three pairs span high to low error correlation:
  a  depression, sleep apnoea (hospital records) vs any sleep disorder  (nested definitions)
  b  PTSD, constipation with vs without laxative purchases              (partly overlapping)
  c  PTSD, EHR vs clinical depression                                   (almost no shared samples)
Axes show pseudovalues centered on each estimate (pseudovalue minus rg), so the three panels share
one scale. Source: results/definition_effect/ldsc/*.delete and rg_common.tsv.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def main():
    vs.apply()
    panels = [
        ("a", "MDD", "G6_SLEEPAPNO", "SLEEP", "MDD", "", "Depression x sleep apnoea"),
        ("b", "PTSD", "K11_CONSTIPATION", "K11_OTHFUNC", "PTSD", "", "PTSD x constipation"),
        ("c", "PTSD", "MDD_EHR", "MDD_Clin", "PTSD", "__control", "PTSD x depression"),
    ]
    label = dict(ENDPOINT_LABEL, MDD_EHR="EHR definition", MDD_Clin="clinical definition")
    fig, axes = plt.subplots(1, 3, figsize=(vs.DOUBLE, 2.75))
    fig.subplots_adjust(left=0.08, right=0.99, top=0.76, bottom=0.2, wspace=0.42)
    lim = 3.6
    for ax, (letter, shared, d1, d2, run, tag, title) in zip(axes, panels):
        p1, r1, s1 = pseudo(shared, d1, run, tag)
        p2, r2, s2 = pseudo(shared, d2, run, tag)
        x, y = p1 - r1, p2 - r2
        n = len(x)
        se_d = np.sqrt(np.var(p1 - p2, ddof=1) / n)
        se_ind = np.sqrt(s1 ** 2 + s2 ** 2)
        corr = np.corrcoef(p1, p2)[0, 1]
        # ribbon: diagonal band holding the middle 90% of per-block differences
        dlo, dhi = np.percentile(y - x, [5, 95])
        t = np.array([-lim, lim])
        ax.fill_between(t, t + dlo, t + dhi, color=vs.TRAIT[shared], alpha=0.12, lw=0, zorder=1)
        ax.plot(t, t, color=vs.MUTED, lw=0.6, zorder=2)
        inside = (np.abs(x) <= lim) & (np.abs(y) <= lim)
        ax.scatter(x[inside], y[inside], s=9, facecolor=vs.TRAIT[shared], edgecolor="white", linewidths=0.3, zorder=3)
        out = int((~inside).sum())
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xticks([-3, 0, 3])
        ax.set_yticks([-3, 0, 3])
        ax.set_xlabel(f"{label[d1]}", fontsize=7.5)
        ax.set_ylabel(f"{label[d2]}", fontsize=7.5)
        ax.set_title(f"{letter}  {title}", loc="left", fontsize=8.5, pad=22)
        ax.text(0.0, 1.035, f"error correlation {corr:.2f}", transform=ax.transAxes, fontsize=7, color=vs.INK_2)
        ax.text(0.03, 0.97, f"SE of difference {se_d:.3f}\nif independent {se_ind:.3f}"
                + (f"\n{out} block{'s' if out > 1 else ''} off scale" if out else ""),
                transform=ax.transAxes, ha="left", va="top", fontsize=7, color=vs.INK_2, zorder=5)
    fig.text(0.08, 0.975, "Each dot is one of 200 genome blocks: its rg pseudovalue under each definition, centered on the estimate.",
             fontsize=7, color=vs.INK_2)
    fig.text(0.08, 0.935, "Tinted band: middle 90% of the per-block differences. Gray line: equal values.", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig11_pseudovalues")


if __name__ == "__main__":
    main()
