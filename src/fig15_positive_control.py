"""Figure 15. Positive control: does the method see a definition effect where one is known?

Depression defined from electronic health records and depression defined by clinical assessment
are two definitions of one condition (rg between them 0.86). For every other trait and FinnGen
endpoint, a dumbbell joins rg with the EHR definition (filled) and with the clinical definition
(hollow). Rows are sorted by the difference within two blocks: psychiatric traits, then medical
endpoints. Filled marks right of hollow marks mean the EHR definition correlates more strongly.
Dark labels mark differences with p < 0.05 from the jackknife test.

Source: results/definition_effect/pairwise_delta.tsv (family "MDD definition (positive control)")
"""
import matplotlib.pyplot as plt
import pandas as pd

import viz_style as vs
from itoju_labels import ENDPOINT_GROUP, ENDPOINT_LABEL, TRAIT_LABEL

ROOT = vs.ROOT
PSY = ["ASD", "SCZ", "BIP", "PTSD"]


def main():
    vs.apply()
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    c = pw[pw.family == "MDD definition (positive control)"].copy()
    c["kind"] = c.shared.isin(PSY).map({True: "Psychiatric traits", False: "Medical endpoints (FinnGen)"})
    def qualified(e):
        # a definition name only makes sense beside its condition ("Epilepsy: generalized")
        if e in TRAIT_LABEL:
            return TRAIT_LABEL[e]
        g, lab = ENDPOINT_GROUP[e], ENDPOINT_LABEL[e]
        if g == "Neighbouring conditions":
            return lab
        if g == "Intellectual disability":
            return "Intellectual disability (F7)"
        if g == "Constipation":
            return lab   # the label already names the condition
        return f"{g}: {lab[0].lower() + lab[1:] if not lab[:2].isupper() else lab}"
    c["label"] = [qualified(e) for e in c.shared]

    blocks = [("Psychiatric traits", c[c.kind == "Psychiatric traits"].sort_values("delta")),
              ("Medical endpoints (FinnGen)", c[c.kind != "Psychiatric traits"].sort_values("delta"))]
    n_rows = sum(len(b) for _, b in blocks) + 2 * len(blocks)
    fig = plt.figure(figsize=(vs.SINGLE, 0.135 * n_rows + 0.9))
    ax = fig.add_axes([0.53, 0.10, 0.45, 0.84])
    col = vs.TRAIT["MDD"]
    y = 0
    for name, b in blocks:
        ax.text(-0.36, y, name, ha="right", va="center", fontsize=7.6, fontweight="bold", color=vs.INK,
                transform=ax.get_yaxis_transform() if False else ax.transData)
        y += 1.2
        for r in b.itertuples():
            sig = r.p < 0.05
            # EHR a little above, clinical a little below, so equal estimates never hide one another
            ax.plot([r.rg2, r.rg1], [y + 0.14, y - 0.14], color=col, lw=1.0, alpha=0.8, zorder=2)
            ax.scatter(r.rg2, y + 0.14, s=20, facecolor="white", edgecolor=col, linewidths=1.0, zorder=3)
            ax.scatter(r.rg1, y - 0.14, s=20, facecolor=col, edgecolor=col, linewidths=1.0, zorder=4)
            ax.text(-0.36, y, r.label + (" *" if sig else ""), ha="right", va="center", fontsize=7,
                    color=vs.INK if sig else vs.INK_2)
            y += 1
        y += 0.8
    ax.axvline(0, color=vs.GRID, lw=0.6, zorder=0)
    ax.set_xlim(-0.3, 1.0)
    ax.set_ylim(y, -1.4)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("rg with depression")
    n_pos = int((c.delta > 0).sum())
    ax.set_title(f"EHR higher in {n_pos} of {len(c)}", loc="left", fontsize=8, color=vs.INK_2, pad=4)
    # key
    ax.scatter(0.34, -0.9, s=20, facecolor=col, edgecolor=col, clip_on=False)
    ax.text(0.38, -0.9, "EHR", va="center", fontsize=7, color=vs.INK_2)
    ax.scatter(0.64, -0.9, s=20, facecolor="white", edgecolor=col, linewidths=1.0, clip_on=False)
    ax.text(0.68, -0.9, "clinical", va="center", fontsize=7, color=vs.INK_2)
    fig.text(0.02, 0.008, "* EHR vs clinical difference p < 0.05 (jackknife)", fontsize=7, color=vs.INK_2)
    vs.save(fig, "fig15_positive_control")


if __name__ == "__main__":
    main()
