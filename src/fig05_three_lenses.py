"""Figure 5. How much do two definitions overlap? It depends on the vocabulary you look through.

For every pair of alternative definitions of the same condition, the overlap (Jaccard index) is
measured four ways, from the codes a registry stores to the clinical meaning behind them:
  ICD-10 codes      the WHO ICD-10 codes each definition includes
  phecodeX          the research groupings those codes fall into
  SNOMED CT         the standard concepts the codes map to in the OMOP vocabularies
  SNOMED + parents  the same concepts widened to their ancestors within two levels
One line per pair, left to right. Pairs with identical profiles are merged and counted. A line
that climbs means two definitions that look disjoint as codes turn out to describe related
clinical ideas; a line that falls means shared codes hide different meanings.

Sources: results/definitions/pairwise_distance.tsv (src/02), results/omop/pairwise_omop.tsv (src/08b).
"""
import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs
from itoju_labels import ENDPOINT_LABEL

ROOT = vs.ROOT
SETS = {
    "Epilepsy": ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"],
    "Sleep apnoea": ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"],
    "Insomnia": ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"],
    "Constipation": ["K11_CONSTIPATION", "K11_OTHFUNC"],
    "ADHD": ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"],
    "Intellectual disability": ["F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE"],
    "Bladder (two conditions)": ["N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD"],
}
LENSES = [("icd10_jaccard", "ICD-10\ncodes"), ("phecode_jaccard", "phecode\nX"),
          ("snomed_jaccard", "SNOMED\nCT"), ("snomed_hier_jaccard", "SNOMED\n+parents")]
SHORT = {"G6_EPLEPSY": "any", "FE": "focal", "FE_STRICT": "focal strict", "FE_MODE": "focal mode",
         "GE": "generalized", "GE_STRICT": "gen. strict", "GE_MODE": "gen. mode"}


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "omop" / "pairwise_omop.tsv", sep="\t")
    rows = []
    for fam, defs in SETS.items():
        sub = pw[pw.def1.isin(defs) & pw.def2.isin(defs)]
        for r in sub.itertuples():
            rows.append({"family": fam, "def1": r.def1, "def2": r.def2,
                         **{k: round(getattr(r, k), 3) for k, _ in LENSES}})
    df = pd.DataFrame(rows)
    keys = [k for k, _ in LENSES]
    groups = df.groupby(["family"] + keys, sort=False)

    fig, ax = plt.subplots(figsize=(vs.SINGLE, 3.6))
    fig.subplots_adjust(left=0.12, right=0.60, top=0.92, bottom=0.13)
    xs = list(range(len(LENSES)))
    highlight = {
        ("Constipation", "K11_CONSTIPATION", "K11_OTHFUNC"): "Constipation vs K59",
        ("ADHD", "F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"): "ADHD F90.0 vs F90",
        ("Epilepsy", "G6_EPLEPSY", "GE"): "Epilepsy any vs gen.",
        ("Insomnia", "F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"): "Insomnia",
        ("Bladder (two conditions)", "N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD"): "Bladder N31 vs N32",
        ("Sleep apnoea", "G6_SLEEPAPNO", "SLEEP"): "Apnoea vs any sleep",
    }
    labels = []
    for (fam, *vals), g in groups:
        members = set(zip(g.def1, g.def2))
        name = next((lab for (f, a, b), lab in highlight.items() if f == fam and ((a, b) in members or (b, a) in members)), None)
        lw, col, z = (1.4, vs.INK, 3) if name else (0.8, vs.MUTED, 2)
        ax.plot(xs, vals, color=col, lw=lw, alpha=0.9 if name else 0.55, zorder=z, solid_capstyle="round")
        ax.scatter(xs, vals, s=14 if name else 8, color=col, zorder=z + 1, edgecolors="white", linewidths=0.5)
        if name:
            labels.append((vals[-1], name + (f" (x{len(g)})" if len(g) > 1 else "")))
    # each label sits at its own line's end height; only nudged apart if two ends are too close
    labels.sort()
    gap, placed = 0.06, []
    for yend, text in labels:
        ly = max(yend, placed[-1] + gap) if placed else yend
        placed.append(ly)
        ax.text(xs[-1] + 0.12, ly, text, fontsize=7, color=vs.INK_2, va="center", clip_on=False)
    ax.text(xs[-1] + 0.12, 1.0, "identical code sets\n(strict and mode variants)", fontsize=7, color=vs.MUTED,
            va="center", clip_on=False)
    ax.set_xticks(xs)
    ax.set_xticklabels([lab for _, lab in LENSES], fontsize=7)
    ax.set_xlim(-0.2, xs[-1] + 0.05)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("overlap between the two definitions (Jaccard)")
    ax.text(-0.2, 1.08, f"{len(df)} pairs of alternative definitions; gray: the rest", transform=ax.transData,
            fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig05_three_lenses")


if __name__ == "__main__":
    main()
