"""Figure: do definitions that differ more give genetic correlations that differ more?

Single column, two panels stacked. Each point is one comparison between two definitions of the same
condition, for one psychiatric trait (135 comparisons), in the colour of that trait; epilepsy's 105 are
hollow. Across is how different the two definitions are; up is how far rg moved.
  a  the pre-specified definition distance: 1 - mean(ICD-10 code Jaccard, data-source Jaccard,
     1 - any rule difference), fixed in the lab notebook before any rg was estimated.
  b  1 - SNOMED CT hierarchy-aware Jaccard (concepts plus ancestors within two levels).
Definition pairs share exact distances, so points are spread sideways by up to 0.012 and fall in clusters;
for each cluster of five or more, a dark bar marks the median |difference| and a pale box the middle half.
Above each panel, a strip places the Spearman correlation for each trait (coloured dots) and for all
comparisons (dark bar) on one scale. The comparisons share traits and definitions, so the points are not
independent and the correlations are descriptive, not tests.

Sources: results/definition_effect/pairwise_delta.tsv, results/omop/pairwise_omop.tsv.
Writes results/definition_effect/distance_vs_delta.tsv.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import itoju_svg as sv

ROOT = sv.ROOT
PANELS = [
    ("distance", "definition distance: codes, sources, rules", "a  Codes, sources and rules", "pre-specified",
     (0.28, 0.88), (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)),
    ("snomed_distance", "SNOMED CT distance (1 - hierarchy Jaccard)", "b  Clinical meaning", "exploratory",
     (-0.04, 1.02), (0, 0.25, 0.5, 0.75, 1.0)),
]
YMAX = 0.5
PLOT_H = 100.0


def main():
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw = pw[pw.family != "MDD definition (positive control)"].copy()
    om = pd.read_csv(ROOT / "results" / "omop" / "pairwise_omop.tsv", sep="\t")
    om = pd.concat([om[["def1", "def2", "snomed_hier_jaccard"]],
                    om.rename(columns={"def1": "def2", "def2": "def1"})[["def1", "def2", "snomed_hier_jaccard"]]])
    pw = pw.merge(om.drop_duplicates(["def1", "def2"]), on=["def1", "def2"], how="left")
    assert pw.snomed_hier_jaccard.notna().all()
    pw["snomed_distance"] = 1 - pw.snomed_hier_jaccard
    pw["abs_delta"] = pw.delta.abs()
    epi = (pw.family == "Epilepsy").values
    rho = {col: spearmanr(pw[col], pw.abs_delta)[0] for col, *_ in PANELS}

    W = sv.SINGLE
    left = 12.0
    ax0, ax1 = 36.0, W - 10
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    title = (["Clinical meaning tracks the shift", "better than code lists do"] if rho["snomed_distance"] > rho["distance"]
             else ["Code lists track the shift", "better than clinical meaning"])
    f.header("itoju  /  distance and drift", title, "|shift in rg| against definition distance")
    f.key_traits(left, 76, ["ASD", "SCZ", "BIP"])
    f.key_traits(left, 88, ["MDD", "PTSD"])
    f.key_row(left, 100, [("dot", sv.INK_2, "sleep, bowel, ADHD"), ("hollow", sv.INK_2, "epilepsy")])
    f.key_row(left, 112, [("line", sv.INK, "cluster median"), ("band", sv.NOISE, "middle half")])

    stats = []
    t0 = 136.0
    jit = np.random.default_rng(7).uniform(-0.012, 0.012, len(pw))
    n_over = 0
    for col, xlabel, head, kind, xlim, xticks in PANELS:
        f.text(left, t0, head, 8.2, sv.INK, family=sv.SERIF)
        f.text(left + sv.text_width(head, 8.2, serif=True) + 6, t0, kind, 7.0, sv.DIM)
        f.text(left, t0 + 10, f"Spearman rho: dots by trait, bar all 135 ({rho[col]:.2f})", 7.0, sv.DIM)
        # rho strip
        S = sv.scale(-0.3, 1.0, ax0, ax1)
        sy = t0 + 24
        f.line(ax0, sy, ax1, sy, sv.GRID, 0.5)
        f.line(S(0), sy - 5, S(0), sy + 5, sv.RULE, 0.5)
        for t in sv.FIVE:
            s = pw[pw.shared == t]
            r_, p_ = spearmanr(s[col], s.abs_delta)
            stats.append({"distance": col, "shared": t, "n": len(s), "spearman_rho": r_, "p_descriptive": p_})
            f.circle(S(r_), sy, 2.4, sv.TRAIT[t], stroke=sv.GROUND, sw=0.5)
        r_all, p_all = spearmanr(pw[col], pw.abs_delta)
        stats.append({"distance": col, "shared": "all", "n": len(pw), "spearman_rho": r_all, "p_descriptive": p_all})
        f.line(S(r_all), sy - 6, S(r_all), sy + 6, sv.INK, 1.8, mark=True)
        for v, s_ in ((0, "0"), (0.5, "0.5"), (1.0, "1")):
            f.text(S(v), sy + 14, s_, 7.0, sv.DIM, anchor="middle")
        f.text(ax0 - 4, sy + 2.4, "rho", 7.0, sv.DIM, anchor="end")

        # scatter
        py0 = sy + 32
        py1 = py0 + PLOT_H
        X = sv.scale(xlim[0], xlim[1], ax0, ax1)
        Y = sv.scale(0, YMAX, py1, py0)
        for v in (0.1, 0.2, 0.3, 0.4, 0.5):
            f.line(ax0, Y(v), ax1, Y(v), sv.GRID, 0.4)
        for v in (0, 0.1, 0.2, 0.3, 0.4, 0.5):
            f.text(ax0 - 4, Y(v) + 2.4, f"{v:g}", 7.0, sv.DIM, anchor="end")
        clusters = np.round(pw[col].values / 0.04)
        boxes = []
        for c_ in np.unique(clusters):
            b = pw[clusters == c_]
            if len(b) < 5:
                continue
            lo, hi = b[col].min() - 0.016, b[col].max() + 0.016
            q1, q3 = b.abs_delta.quantile(0.25), b.abs_delta.quantile(0.75)
            f.rect(X(lo), Y(q3), X(hi) - X(lo), Y(q1) - Y(q3), sv.NOISE)
            boxes.append((X(lo), X(hi), Y(b.abs_delta.median())))
        for i, r in enumerate(pw.itertuples()):
            x, yv = X(getattr(r, col) + jit[i]), r.abs_delta
            if yv > YMAX:
                n_over += 1
                yv = YMAX
            if epi[i]:
                f.circle(x, Y(yv), 1.5, "none", stroke=sv.TRAIT[r.shared], sw=0.7)
            else:
                f.circle(x, Y(yv), 1.8, sv.TRAIT[r.shared], stroke=sv.GROUND, sw=0.3)
        for xa, xb, ym in boxes:
            f.line(xa, ym, xb, ym, sv.INK, 1.8, mark=True)
        f.line(ax0, py1, ax1, py1, sv.RULE, 0.6)
        for v in xticks:
            f.line(X(v), py1, X(v), py1 + 2.5, sv.RULE, 0.5)
            f.text(X(v), py1 + 10.5, f"{v:g}", 7.0, sv.DIM, anchor="middle")
        f.text((ax0 + ax1) / 2, py1 + 20, xlabel, 7.0, sv.DIM, anchor="middle")
        f.text(left, py0 - 6, "|difference in rg|", 7.0, sv.DIM)
        t0 = py1 + 40
    assert n_over == 0, f"{n_over} comparisons above the axis"

    st = pd.DataFrame(stats)
    st.to_csv(ROOT / "results" / "definition_effect" / "distance_vs_delta.tsv", sep="\t", index=False)
    print(st.round(3).to_string(index=False))
    by = lambda col: st[(st.distance == col) & (st.shared != "all")].spearman_rho
    lines = [f"Codes, sources, rules: rho {rho['distance']:.2f};",
             f"SNOMED CT meaning: rho {rho['snomed_distance']:.2f}.",
             f"By trait, {by('distance').min():.2f} to {by('distance').max():.2f} against "
             f"{by('snomed_distance').min():.2f} to {by('snomed_distance').max():.2f}.",
             "Points share traits and definitions:",
             "these correlations describe, not test."]
    end = f.reading(left, t0 - 6, W - left - 10, lines, strong=(0, 1))
    f.h = end + 20
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.source("pairwise_delta.tsv, pairwise_omop.tsv")
    f.save("fig14_distance_drift")


if __name__ == "__main__":
    main()
