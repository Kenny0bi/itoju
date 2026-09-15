"""Figure 20. What buys the power to see a definition effect: shared cases.

The smallest detectable difference between two genetic correlations is
    MDD = (z_0.975 + z_0.80) * sqrt(se1^2 + se2^2 - 2 * rho * se1 * se2),
where rho is the correlation between the two estimates' sampling errors. The shaded map is that
formula for two estimates with the same SE (vertical axis) and error correlation rho (horizontal
axis), drawn like a topographic map in warm grays, because it is the reference landscape: each contour
line is one detectable difference. Two definitions that share most of their cases sit far right, where
the contours plunge.

Each dot is one comparison between definitions of the same condition (135), in the color of its
psychiatric trait, placed at its jackknife error correlation and the root mean square of its two SEs;
turmeric rings are the 31 EHR vs clinical depression comparisons. With unequal SEs the map is an
approximation; the script reports how many dots fall in the band of their exact value. Four
comparisons are numbered and named below the plot, so no label covers the data.

Sources: results/definition_effect/pairwise_delta.tsv and rg_common.tsv.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

import viz_style as vs
from itoju_labels import TRAIT_SHORT

ROOT = vs.ROOT
K = norm.ppf(0.975) + norm.ppf(0.80)
LEVELS = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 10.0]
BANDS = ["#f7f5f0", "#ece8df", "#ddd7cb", "#c9c1b2", "#aea596", "#8b8376", "#686156"]
YLIM = (0.012, 0.7)
FIVE = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
EXAMPLES = [  # shared trait, def1, def2, name, offset of the number (points)
    ("MDD", "G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "sleep apnoea, hospital vs + primary care (depression)", (-9, 4)),
    ("PTSD", "K11_CONSTIPATION", "K11_OTHFUNC", "constipation or laxatives vs all of K59 (PTSD)", (0, 9)),
    ("MDD", "FE", "GE", "focal vs generalized epilepsy (depression)", (0, 9)),
    ("PTSD", "MDD_EHR", "MDD_Clin", "EHR vs clinical depression (PTSD)", (0, -9)),
]


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    rc = pd.read_csv(ROOT / "results" / "definition_effect" / "rg_common.tsv", sep="\t")
    se = {}
    for r in rc.itertuples():
        se[(r.trait1, r.trait2)] = se[(r.trait2, r.trait1)] = r.se
    pw["se1"] = [se[(r.shared, r.def1)] for r in pw.itertuples()]
    pw["se2"] = [se[(r.shared, r.def2)] for r in pw.itertuples()]
    pw["se_rms"] = np.sqrt((pw.se1 ** 2 + pw.se2 ** 2) / 2)
    exact = K * np.sqrt(pw.se1 ** 2 + pw.se2 ** 2 - 2 * pw.error_corr * pw.se1 * pw.se2)
    print(f"  exact formula vs stored min_detectable_delta: max abs diff {np.max(np.abs(exact - pw.min_detectable_delta)):.4f}")
    band_of = lambda v: np.searchsorted(LEVELS, v) - 1
    same_band = np.mean(band_of(exact.values) == band_of(K * pw.se_rms.values * np.sqrt(2 - 2 * pw.error_corr.values)))
    print(f"  dots whose map band matches their exact band: {same_band:.1%}")

    fig, ax = plt.subplots(figsize=(vs.SINGLE, 4.3))
    fig.subplots_adjust(left=0.165, right=0.975, top=0.745, bottom=0.265)
    rho = np.linspace(-0.2, 0.999, 400)
    s = np.logspace(np.log10(YLIM[0]), np.log10(YLIM[1]), 400)
    R, S = np.meshgrid(rho, s)
    M = K * S * np.sqrt(2 - 2 * R)
    ax.contourf(R, S, M, levels=LEVELS, colors=BANDS, zorder=0)
    ax.contour(R, S, M, levels=LEVELS[1:-1], colors="white", linewidths=0.7, zorder=1)
    # each contour labelled where it enters the view; ink on the light bands, white on the two darkest
    x_left, y_floor = -0.19, YLIM[0] * 1.12
    for lv in LEVELS[1:-1]:
        y_lab = lv / (K * np.sqrt(2 - 2 * x_left))
        ink = vs.INK if lv <= 0.2 else "white"
        if YLIM[0] * 1.15 < y_lab < YLIM[1] * 0.85:
            ax.text(x_left, y_lab * 1.04, f"{lv:g}", ha="left", va="bottom", fontsize=7, color=ink, zorder=2)
        else:
            r_lab = 1 - (lv / (K * y_floor)) ** 2 / 2
            if -0.15 < r_lab < 0.97:
                ax.text(r_lab + 0.015, y_floor, f"{lv:g}", ha="left", va="bottom", fontsize=7, color=ink, zorder=2)

    ctrl = pw.family == "MDD definition (positive control)"
    main_ = pw[~ctrl]
    ax.scatter(main_.error_corr, main_.se_rms, s=13, c=[vs.TRAIT[t] for t in main_.shared], edgecolor="white",
               linewidths=0.4, zorder=3)
    ax.scatter(pw.error_corr[ctrl], pw.se_rms[ctrl], s=12, facecolor="white", edgecolor=vs.TRAIT["MDD"],
               linewidths=0.9, zorder=3)
    lines = []
    for k, (shared, d1, d2, name, off) in enumerate(EXAMPLES, start=1):
        r = pw[(pw.shared == shared) & (pw.def1 == d1) & (pw.def2 == d2)].iloc[0]
        print(f"  example {k}: rho {r.error_corr:.3f}, se_rms {r.se_rms:.4f}, MDD {r.min_detectable_delta:.3f}")
        ax.scatter(r.error_corr, r.se_rms, s=34, facecolor="none", edgecolor=vs.INK, linewidths=0.9, zorder=4)
        ax.annotate(str(k), (r.error_corr, r.se_rms), xytext=off, textcoords="offset points", ha="center",
                    va="center", fontsize=7, color=vs.INK, zorder=5,
                    bbox=dict(boxstyle="circle,pad=0.18", facecolor="white", edgecolor=vs.INK, lw=0.7))
        lines.append(f"{k}  {name}: detects {r.min_detectable_delta:.2f}")
    ax.set_yscale("log")
    ax.set_ylim(*YLIM)
    ax.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5])
    ax.set_yticklabels(["0.02", "0.05", "0.1", "0.2", "0.5"])
    ax.minorticks_off()
    ax.set_xlim(-0.2, 1.0)
    ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("how much the two estimates move together\n(error correlation; high = shared cases)", fontsize=7.5)
    ax.set_ylabel("SE of each estimate", fontsize=7.5)
    fig.text(0.02, 0.985, "Gray map and contours: the smallest difference in rg\nthe test would catch 4 times in 5.",
             fontsize=7.5, color=vs.INK, va="top", linespacing=1.3)
    fig.text(0.02, 0.9, "Dot: one comparison of two definitions (135), in its trait's color.", fontsize=7,
             color=vs.INK_2, va="top")
    x = 0.02
    for t in FIVE:
        fig.add_artist(plt.Line2D([x + 0.012], [0.835], marker="o", markersize=4.2, color=vs.TRAIT[t], ls="",
                                  transform=fig.transFigure))
        fig.text(x + 0.03, 0.835, TRAIT_SHORT[t], fontsize=7, color=vs.INK_2, va="center")
        x += 0.05 + 0.024 * len(TRAIT_SHORT[t])
    fig.add_artist(plt.Line2D([0.032], [0.8], marker="o", markersize=4.2, markerfacecolor="white",
                              markeredgecolor=vs.TRAIT["MDD"], ls="", transform=fig.transFigure))
    fig.text(0.05, 0.8, "ring: EHR vs clinical depression (31)", fontsize=7, color=vs.INK_2, va="center")
    fig.text(0.02, 0.105, "\n".join(lines), fontsize=7, color=vs.INK_2, va="top", linespacing=1.35)
    vs.save(fig, "fig20_power")


if __name__ == "__main__":
    main()
