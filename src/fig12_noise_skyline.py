"""Figure: every definition change in the study, inside its own noise.

One vertical mark per comparison (166): 135 between definitions of the same condition for five psychiatric
GWAS, and 31 between depression defined from health records and by clinical assessment. The gray sleeve
spans plus or minus 1.96 SE of the difference from the shared-block jackknife; the stem is the shift in
genetic correlation (second definition minus first; for depression, EHR minus clinical). A stem that
escapes its sleeve (p < 0.05) is drawn thick with an end dot in its trait's colour; one that stays inside
is thin. Families get different widths per comparison so the small ones stay readable beside epilepsy's
105. Sleep apnoea's shifts are small and precise, so its panel has its own labelled scale.

Source: results/definition_effect/pairwise_delta.tsv. Every count in the figure is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv
from itoju_labels import ENDPOINT_ORDER

ROOT = sv.ROOT
FAMILIES = [  # family, title, (unused), y-axis half range
    ("Sleep apnoea", "Sleep apnoea", 7.2, 0.08),
    ("Insomnia", "Insomnia", 7.2, 0.5),
    ("Constipation", "Constipation", 7.2, 0.5),
    ("ADHD", "ADHD", 7.2, 0.5),
    ("Epilepsy", "Epilepsy", 1.34, 0.5),
    ("MDD definition (positive control)", "EHR vs clinical", 1.34, 0.5),
]
CODES = {"G6_SLEEPAPNO": "A", "G6_SLEEPAPNO_INCLAVO": "B", "SLEEP": "C"}
PAIR_GAP = 0.8
TOP, HEIGHT = 130.0, 206.0
GAP, GAP_SCALE, AXIS_W = 9.0, 26.0, 28.0


def ordered(sub, fam):
    rank = {s: i for i, s in enumerate(sv.FIVE + ENDPOINT_ORDER)}
    if fam.startswith("MDD"):
        return [list(sub.sort_values("shared", key=lambda s: s.map(rank)).itertuples())]
    return [list(g.sort_values("shared", key=lambda s: s.map(rank)).itertuples())
            for _, g in sub.groupby(["def1", "def2"], sort=False)]


def main():
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    pw["shift"] = -pw.delta                                  # delta = rg(def1) - rg(def2)
    ctrl = pw.family == "MDD definition (positive control)"
    assert (pw[ctrl].def1 == "MDD_EHR").all()
    pw.loc[ctrl, "shift"] = pw.loc[ctrl, "delta"]            # EHR minus clinical, so up means EHR higher
    pw["hw"] = 1.96 * pw.se_delta
    pw["clear"] = pw.p < 0.05

    layouts = []
    for fam, title, unit, half in FAMILIES:
        groups = ordered(pw[pw.family == fam], fam)
        xs, x = [], 0.0
        for g in groups:
            for r in g:
                xs.append((x, r))
                x += 1
            x += PAIR_GAP
        layouts.append(dict(fam=fam, title=title, half=half, groups=groups, xs=xs, span=x - PAIR_GAP))
    W = sv.DOUBLE
    left = 14.0
    ROWS = [(["Sleep apnoea", "Insomnia", "Constipation", "ADHD"], 138.0), (["Epilepsy", "MDD definition (positive control)"], 350.0)]
    HP = 132.0
    H = 356.0 + HP + 100
    f = sv.Figure(W, H)
    f.header("itoju  /  every comparison",
             "Every definition change in the study, inside its own noise",
             f"{int((~ctrl).sum())} changes of definition for five psychiatric GWAS, and {int(ctrl.sum())} between EHR and clinical depression")
    f.key_traits(left, 72)
    f.key_row(left, 86, [("band", sv.NOISE, "noise, 1.96 SE of the difference"),
                         ("line", sv.INK_2, "escapes its noise (p < 0.05): thick, with a dot"),
                         ("line", sv.RULE, "inside it: thin")])
    by = {l["fam"]: l for l in layouts}
    counts = {}
    for names, top in ROWS:
        row = [by[n] for n in names]
        gaps = [0.0] + [GAP_SCALE if row[k]["half"] != row[k - 1]["half"] else GAP for k in range(1, len(row))]
        raw = [l["span"] + 1 for l in row]
        k_w = (W - left - 10 - AXIS_W - sum(gaps)) / sum(raw)
        y0 = top + HP / 2
        x_cursor = left + AXIS_W
        for k, l in enumerate(row):
            x_cursor += gaps[k]
            half = l["half"]
            wpt = raw[k] * k_w
            u = wpt / (l["span"] + 1)
            X = lambda v, x0=x_cursor, u=u: x0 + (v + 0.5) * u
            Y = lambda v, half=half, y0=y0: y0 - v / half * HP / 2
            if k == 0 or gaps[k] == GAP_SCALE:
                ticks = (-0.05, 0, 0.05) if half < 0.1 else (-0.4, -0.2, 0, 0.2, 0.4)
                axx = x_cursor - 3
                f.line(axx, top, axx, top + HP, sv.RULE, 0.5)
                for v in ticks:
                    f.line(axx - 2.5, Y(v), axx, Y(v), sv.RULE, 0.5)
                    f.text(axx - 4, Y(v) + 2.4, "0" if v == 0 else sv.fmt(v, 2 if half < 0.1 else 1, sign=True), 7.0,
                           sv.DIM, anchor="end")
            f.line(x_cursor - 1, y0, x_cursor + wpt + 1, y0, sv.INK_2, 0.5)
            n_clear = 0
            for pos, r in l["xs"]:
                hwc = min(r.hw, half)
                sw = max(0.4, 0.38 * u)
                f.rect(X(pos) - sw, Y(hwc), 2 * sw, Y(-hwc) - Y(hwc), sv.NOISE)
                col = sv.TRAIT["MDD" if l["fam"].startswith("MDD") else r.shared]
                s = max(min(r.shift, half), -half)
                if r.clear:
                    n_clear += 1
                    f.line(X(pos), y0, X(pos), Y(s), col, max(0.9, min(2.4, 0.3 * u)), mark=True)
                    f.circle(X(pos), Y(s), max(1.1, min(2.6, 0.34 * u)), col)
                else:
                    f.line(X(pos), y0, X(pos), Y(s), col, max(0.35, min(1.0, 0.14 * u)), mark=True)
            counts[l["fam"]] = (n_clear, len(l["xs"]))
            f.text(x_cursor, top - 20, l["title"], 8.2, sv.INK, family=sv.SERIF)
            f.text(x_cursor, top - 10, f"{n_clear} of {len(l['xs'])} escape" + ("   own scale" if half < 0.1 else ""), 7.0, sv.DIM)
            yb = top + HP + 12
            fam = l["fam"]
            if fam == "Sleep apnoea":
                for g in l["groups"]:
                    gx = np.mean([X(p) for p, r in l["xs"] if (r.def1, r.def2) == (g[0].def1, g[0].def2)])
                    f.text(gx, yb, f"{CODES[g[0].def1]} to {CODES[g[0].def2]}", 7.0, sv.INK_2, anchor="middle")
                notes = ["A hospital records", "B + primary care", "C any sleep disorder"]
                yb += 11
            elif fam == "Insomnia":
                notes = ["F51.0 to", "all of F51"]
            elif fam == "Constipation":
                notes = ["K59.0 or", "laxatives to", "all of K59"]
            elif fam == "ADHD":
                notes = ["F90.0 to", "all of F90"]
            elif fam == "Epilepsy":
                sub = pw[pw.family == fam]
                kind = lambda d: "F" if d.startswith("FE") else ("G" if d.startswith("GE") else "A")
                k1, k2 = sub.def1.map(kind), sub.def2.map(kind)
                fg = (k1 != k2) & (k1 != "A") & (k2 != "A")
                same = k1 == k2
                anyv = ~fg & ~same
                counts["epi_fg"] = (int(sub[fg].clear.sum()), int(fg.sum()))
                counts["epi_any"] = (int(sub[anyv].clear.sum()), int(anyv.sum()))
                counts["epi_same"] = (int(sub[same].clear.sum()), int(same.sum()))
                notes = [f"21 pairs of the 7 definitions, five traits each: focal vs generalized {counts['epi_fg'][0]} of "
                         f"{counts['epi_fg'][1]} escape,",
                         f"any vs a subtype {counts['epi_any'][0]} of {counts['epi_any'][1]}, strict or mode vs plain "
                         f"{counts['epi_same'][0]} of {counts['epi_same'][1]}"]
            else:
                notes = ["depression: EHR minus", "clinical, 31 traits"]
            for i, s in enumerate(notes):
                f.text(x_cursor, yb + i * 9, s, 7.0, sv.DIM, max_w=wpt + gaps[min(k + 1, len(gaps) - 1)] - 3)
            x_cursor += wpt

    print("  escapes:", counts)
    sa, ins, con, adhd, epi, cc = (counts[fam] for fam, *_ in FAMILIES)
    lines = [
        f"Widening the scope: sleep apnoea {sa[0]} of {sa[1]} escape, insomnia {ins[0]} of {ins[1]}, "
        f"constipation {con[0]} of {con[1]}; ADHD {adhd[0]} of {adhd[1]}.",
        f"Epilepsy: focal against generalized {counts['epi_fg'][0]} of {counts['epi_fg'][1]}; strict or mode against "
        f"plain {counts['epi_same'][0]} of {counts['epi_same'][1]}.",
        f"EHR against clinical depression: {cc[0]} of {cc[1]}; the two GWAS share almost no people, so sleeves stay wide.",
    ]
    f.reading(left, 350.0 + HP + 34, W - left - 10, lines, strong=(0,))
    f.source("shared-block jackknife on the common SNP set; stems are rg(second) minus rg(first); sleeves cut at the axis ends")
    f.save("fig12_noise_skyline")


if __name__ == "__main__":
    main()
