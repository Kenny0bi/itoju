"""Figure 15. Positive control: depression defined from health records casts a wider genetic net.

Depression defined from electronic health records and depression defined by clinical assessment are
two definitions of one condition (rg between them 0.86). The known result (Cai et al. 2020) is that
the looser definition carries less specific genetics, so it should overlap more with everything else.

Left: one spoke per trait or FinnGen endpoint (31), grouped by body system. Distance from the center
is rg with depression. The filled shape joins rg with the EHR definition; the dark outline joins rg
with the clinical definition. Where the fill reaches past the outline, the EHR definition correlates
more strongly. Spokes whose difference passes the jackknife test (p < 0.05) are drawn in color.
Right: the same 31 differences (EHR minus clinical), ranked, each inside its noise sleeve (plus or
minus 1.96 SE of the difference); a colored stem escapes its sleeve.

Source: results/definition_effect/pairwise_delta.tsv (family "MDD definition (positive control)").
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Polygon

import viz_style as vs

ROOT = vs.ROOT
SLEEVE = "#e6e1d6"
GRAY = "#b9b5ac"
TINT = "#f8e7bb"
COL = vs.TRAIT["MDD"]
RMIN, RMAX, HOLE = -0.2, 1.0, 0.12
XCAP = 0.6
# spokes in order, grouped by body system; a gap separates the groups
GROUPS = [
    [("ASD", "Autism"), ("SCZ", "Schizophrenia"), ("BIP", "Bipolar"), ("PTSD", "PTSD")],
    [("F5_ADHD", "ADHD (F90.0)"), ("KRA_PSY_HYPERKIN_EXMORE", "ADHD (F90)"),
     ("KRA_PSY_MENTALRET_EXMORE", "Intellectual disability"), ("KRA_PSY_DEVWIDE_EXMORE", "Autism spectrum (F84)"),
     ("KRA_PSY_AUTISM_EXMORE", "Autism (F84.0, F84.5)")],
    [("G6_EPLEPSY", "Epilepsy, any"), ("FE", "Focal epilepsy"), ("FE_STRICT", "Focal, strict"),
     ("FE_MODE", "Focal, mode"), ("GE", "Generalized epilepsy"), ("GE_STRICT", "Generalized, strict"),
     ("GE_MODE", "Generalized, mode"), ("G6_STATUSEPI", "Status epilepticus")],
    [("G6_SLEEPAPNO", "Apnoea, hospital"), ("G6_SLEEPAPNO_INCLAVO", "Apnoea + primary care"),
     ("SLEEP", "Any sleep disorder"), ("F5_INSOMNIA", "Insomnia"), ("KRA_PSY_SLEEP_NONORG_EXMORE", "Nonorganic sleep"),
     ("G6_SLEEPDISOTH", "Other sleep disorders"), ("F5_SLEEP_NOS", "Sleep, unspecified")],
    [("K11_CONSTIPATION", "Constipation/laxatives"), ("K11_OTHFUNC", "Functional bowel (K59)"),
     ("K11_IBS", "Irritable bowel"), ("K11_FUNCDYSP", "Functional dyspepsia"), ("K11_REFLUX", "Reflux")],
    [("N14_NEUROMUSCDYSBLADD", "Neurogenic bladder"), ("N14_OTHBLADD", "Other bladder")],
]


def radius(v):
    return HOLE + (np.clip(v, RMIN, RMAX) - RMIN) / (RMAX - RMIN) * (1 - HOLE)


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    c = pw[pw.family == "MDD definition (positive control)"].copy()
    ehr_first = c.def1 == "MDD_EHR"
    c["ehr"] = np.where(ehr_first, c.rg1, c.rg2)
    c["clin"] = np.where(ehr_first, c.rg2, c.rg1)
    c["diff"] = c.ehr - c.clin
    c["hw"] = 1.96 * c.se_delta
    c["clear"] = c["diff"].abs() > c.hw
    c = c.set_index("shared")
    spokes = [s for g in GROUPS for s in g]
    assert sorted(k for k, _ in spokes) == sorted(c.index), "every control comparison needs a spoke"
    label = dict(spokes)

    fig = plt.figure(figsize=(vs.DOUBLE, 4.45))
    ax = fig.add_axes([0.01, 0.0, 0.60, 0.875])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.8, 1.8)

    slots = sum(len(g) for g in GROUPS) + len(GROUPS)
    theta, pos = {}, 0
    for g in GROUPS:
        for key, _ in g:
            theta[key] = np.deg2rad(90 - 360 * (pos + 0.5) / slots)
            pos += 1
        pos += 1
    gap_angle = np.deg2rad(90 - 360 * (slots - 0.5) / slots)

    for v, lw in ((0.0, 0.8), (0.5, 0.5), (1.0, 0.5)):
        ax.add_patch(Circle((0, 0), radius(v), fill=False, edgecolor=vs.GRID if v else "#cfccc3", lw=lw, zorder=0))
        # ring labels sit in the gap between the last and first group, above the filled shape
        ax.text(radius(v) * np.cos(gap_angle), radius(v) * np.sin(gap_angle), f"{v:g}", ha="center", va="center",
                fontsize=7, color=vs.INK_2, zorder=7,
                bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none"))
    for key, _ in spokes:
        t = theta[key]
        ax.plot([HOLE * np.cos(t), np.cos(t)], [HOLE * np.sin(t), np.sin(t)], color=vs.GRID, lw=0.4, zorder=0)

    def ring(col):
        return np.array([[radius(c.loc[k, col]) * np.cos(theta[k]), radius(c.loc[k, col]) * np.sin(theta[k])]
                         for k, _ in spokes])
    e, cl = ring("ehr"), ring("clin")
    ax.add_patch(Polygon(e, closed=True, facecolor=TINT, edgecolor=COL, lw=1.1, zorder=2))
    ax.add_patch(Polygon(cl, closed=True, fill=False, edgecolor=vs.INK_2, lw=0.9, zorder=3))
    ax.scatter(cl[:, 0], cl[:, 1], s=7, facecolor="white", edgecolor=vs.INK_2, linewidths=0.6, zorder=4)
    for (key, name), pe, pc in zip(spokes, e, cl):
        t = theta[key]
        clear = bool(c.loc[key, "clear"])
        if clear:
            ax.plot([pc[0], pe[0]], [pc[1], pe[1]], color=COL, lw=2.4, solid_capstyle="butt", zorder=5)
            ax.scatter(*pe, s=16, color=COL, linewidths=0, zorder=6)
        deg = np.rad2deg(t)
        right = np.cos(t) >= -1e-9
        r_lab = 1.04
        ax.text(r_lab * np.cos(t), r_lab * np.sin(t), name, rotation=deg if right else deg + 180,
                rotation_mode="anchor", ha="left" if right else "right", va="center", fontsize=7,
                color=vs.INK if clear else vs.INK_2, fontweight="bold" if clear else "normal")

    # key for the net
    fig.text(0.02, 0.965, "Filled shape: rg with depression defined from health records", fontsize=7.5, color=vs.INK_2)
    fig.text(0.02, 0.935, "Dark outline: rg with depression defined by clinical assessment", fontsize=7.5, color=vs.INK_2)
    fig.text(0.02, 0.905, "Rings: rg 0, 0.5, 1.  Colored spoke, bold name: difference p < 0.05", fontsize=7.5,
             color=vs.INK_2)

    # ranked differences, each in its sleeve
    order = c.sort_values("diff", ascending=False)
    n = len(order)
    b = fig.add_axes([0.785, 0.10, 0.205, 0.775])
    ys = np.arange(n)
    hw = np.minimum(order.hw.values, XCAP)
    b.hlines(ys, -hw, hw, colors=SLEEVE, linewidth=4.2, zorder=1)
    b.axvline(0, color=vs.INK_2, lw=0.5, zorder=2)
    cols = [COL if k else GRAY for k in order.clear]
    b.hlines(ys, 0, order["diff"].values, colors=cols, linewidth=1.1, zorder=3)
    k = order.clear.values
    b.scatter(order["diff"].values[k], ys[k], s=12, color=COL, linewidths=0, zorder=4)
    for y, (key, r) in zip(ys, order.iterrows()):
        b.text(-XCAP - 0.04, y, label[key], ha="right", va="center", fontsize=7,
               color=vs.INK if r.clear else vs.INK_2, fontweight="bold" if r.clear else "normal")
    b.set_ylim(n - 0.4, -0.6)
    b.set_xlim(-XCAP - 0.02, XCAP + 0.02)
    b.set_yticks([])
    b.spines["left"].set_visible(False)
    b.set_xticks([-0.5, 0, 0.5])
    b.set_xticklabels(["-0.5", "0", "+0.5"])
    b.set_xlabel("EHR minus clinical, rg", fontsize=7.5)
    n_pos, n_clear = int((c["diff"] > 0).sum()), int(c.clear.sum())
    fig.text(0.785 - 0.155, 0.955, f"EHR higher in {n_pos} of {n}; {n_clear} escape their sleeve", fontsize=8,
             color=vs.INK)
    fig.text(0.785 - 0.155, 0.925, "Sleeve: 1.96 SE of the difference, cut at 0.6", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig15_positive_control")


if __name__ == "__main__":
    main()
