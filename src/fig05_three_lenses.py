"""Figure 5. How much do two definitions overlap? It depends on the vocabulary you look through.

For every pair of alternative definitions of the same condition, overlap is measured four ways, from
the codes a registry stores to the clinical meaning behind them:
  ICD-10 codes      the WHO ICD-10 codes each definition includes
  phecodeX          the research groupings those codes fall into
  SNOMED CT         the standard concepts the codes map to in the OMOP vocabularies
  SNOMED + parents  the same concepts widened to their ancestors within two levels
Each overlap is drawn as an eclipse: two equal discs placed so that the dark shared area, divided by
the area the two discs cover together, equals the Jaccard index. Apart means nothing shared; one disc
means identical. Pairs with identical profiles across all four lenses share a row.

Sources: results/omop/pairwise_omop.tsv (src/08b).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Polygon

import viz_style as vs

ROOT = vs.ROOT
SETS = {
    "Epilepsy": ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"],
    "Sleep apnoea": ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"],
    "Insomnia": ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"],
    "Constipation": ["K11_CONSTIPATION", "K11_OTHFUNC"],
    "ADHD": ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"],
    "Intellectual disability": ["F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE"],
    "Bladder": ["N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD"],
}
SHORT = {"G6_EPLEPSY": "any", "FE": "focal", "FE_STRICT": "focal strict", "FE_MODE": "focal mode",
         "GE": "generalized", "GE_STRICT": "generalized strict", "GE_MODE": "generalized mode",
         "G6_SLEEPAPNO": "hospital", "G6_SLEEPAPNO_INCLAVO": "+ primary care", "SLEEP": "any sleep disorder",
         "F5_INSOMNIA": "F51.0 or G47.0", "KRA_PSY_SLEEP_NONORG_EXMORE": "all of F51",
         "K11_CONSTIPATION": "K59.0 or laxatives", "K11_OTHFUNC": "all of K59", "F5_ADHD": "F90.0",
         "KRA_PSY_HYPERKIN_EXMORE": "all of F90", "F5_MILDRET": "mild (F70)", "KRA_PSY_MENTALRET_EXMORE": "any (F7)",
         "N14_NEUROMUSCDYSBLADD": "neurogenic (N31)", "N14_OTHBLADD": "other (N32)"}
LENSES = [("icd10_jaccard", "ICD-10\ncodes"), ("phecode_jaccard", "phecodeX\ngroups"),
          ("snomed_jaccard", "SNOMED CT\nconcepts"), ("snomed_hier_jaccard", "SNOMED CT\n+ parents")]
R = 0.38
DISC = "#ece9e1"
FIRST = "#c9761a"    # first definition of the pair; validated with SECOND
SECOND = "#2f6fc4"
SHARED = "#1d1d2b"


def lens_area(d, r=R):
    if d >= 2 * r:
        return 0.0
    return 2 * r * r * np.arccos(d / (2 * r)) - (d / 2) * np.sqrt(4 * r * r - d * d)


def separation(j, r=R):
    """Distance between centers of two equal discs whose shared area over covered area equals j."""
    if j <= 0:
        return 2.25 * r
    if j >= 1:
        return 0.0
    lo, hi = 0.0, 2 * r
    for _ in range(60):
        mid = (lo + hi) / 2
        a = lens_area(mid, r)
        if a / (2 * np.pi * r * r - a) > j:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def eclipse(ax, cx, cy, j):
    d = separation(j)
    ca, cb = cx - d / 2, cx + d / 2
    for c, col in ((ca, FIRST), (cb, SECOND)):
        ax.add_patch(Circle((c, cy), R, facecolor=col, edgecolor="white", lw=0.5, zorder=2))
    if j >= 1:
        ax.add_patch(Circle((cx, cy), R, facecolor=SHARED, edgecolor=SHARED, lw=0.6, zorder=3))
    elif d < 2 * R:
        th = np.arccos((d / 2) / R)
        t1 = np.linspace(-th, th, 40)
        t2 = np.linspace(np.pi - th, np.pi + th, 40)
        pts = np.vstack([np.c_[ca + R * np.cos(t1), cy + R * np.sin(t1)],
                         np.c_[cb + R * np.cos(t2), cy + R * np.sin(t2)]])
        ax.add_patch(Polygon(pts, closed=True, facecolor=SHARED, edgecolor="none", zorder=3))
    # check the drawn geometry against the value it stands for
    if 0 < j < 1:
        a = lens_area(d)
        assert abs(a / (2 * np.pi * R * R - a) - j) < 1e-3


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "omop" / "pairwise_omop.tsv", sep="\t")
    keys = [k for k, _ in LENSES]
    rows = []
    for fam, defs in SETS.items():
        sub = pw[pw.def1.isin(defs) & pw.def2.isin(defs)].copy()
        sub[keys] = sub[keys].round(3)
        for vals, g in sub.groupby(keys, sort=False):
            first = g.iloc[0]
            label = f"{SHORT[first.def1]} vs {SHORT[first.def2]}"
            if len(g) > 1:
                label += f"  (+{len(g) - 1} pairs alike)"
            rows.append((fam, label, list(vals)))
    n = len(rows) + len(SETS) + 3
    fig = plt.figure(figsize=(vs.DOUBLE, 0.172 * n + 0.6))
    ax = fig.add_axes([0.30, 0.02, 0.69, 0.98 - 0.45 / (0.172 * n + 0.6)])
    # equal aspect makes one unit the height of a row, so columns sit several units apart
    colw = 5.0
    ax.set_xlim(-0.8, len(LENSES) * colw - 1.5)
    y, last = 0.0, None
    for fam, label, vals in rows:
        if fam != last:
            ax.text(-0.62, y, fam, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK)
            y += 1.0
            last = fam
        ax.text(-0.62, y, label, ha="right", va="center", fontsize=7, color=vs.INK_2)
        for k, v in enumerate(vals):
            eclipse(ax, k * colw + 0.6, y, v)
            ax.text(k * colw + 1.95, y, f"{v:.2f}", ha="left", va="center", fontsize=7, color=vs.INK_2)
        y += 1.0
    ax.set_ylim(y - 0.4, -3.7)
    ax.set_aspect("equal")
    ax.axis("off")
    for k, (_, name) in enumerate(LENSES):
        ax.text(k * colw + 1.1, -1.05, name, ha="center", va="center", fontsize=7.5, color=vs.INK)
    # key: three sample eclipses, drawn by the same function
    ax.text(-0.62, -2.9, "How to read an eclipse", ha="right", va="center", fontsize=7.5, color=vs.INK)
    # spaced by label length, so no disc sits on the label before it
    for x, v, lab in ((0.6, 0.0, "0: nothing shared"), (7.3, 1 / 3, "0.33: a third"), (13.4, 1.0, "1: identical")):
        eclipse(ax, x, -2.9, v)
        ax.text(x + 1.35, -2.9, lab, ha="left", va="center", fontsize=7, color=vs.INK_2)
    fig.text(0.02, 0.99, "Ochre disc: first definition.  Blue disc: second.  Dark area: what they share, as a share of everything either includes (Jaccard).",
             fontsize=7, color=vs.INK_2, va="top")
    vs.save(fig, "fig05_three_lenses")


if __name__ == "__main__":
    main()
