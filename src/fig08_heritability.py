"""Figure: does a wider net lower heritability? SNP heritability against the number of cases.

A definition that casts a wider net gains cases, and with them precision, but it can also pull in people
whose condition is less genetic. One panel per condition; each definition sits at its number of cases (log
scale, each panel its own range) and its liability-scale SNP heritability. The dot is the estimate, the
bar through it the 50% interval, the pale sleeve the 95% interval (cut at zero). The condition's first
definition is ochre, the alternatives blue, the same roles as the other definition figures. A hollow dot
on the floor is an estimate below zero. All panels share the heritability axis.

The liability conversion uses the sample prevalence as the population prevalence, a stated approximation
(FinnGen is a biobank, not a random population sample), so heights compare within a condition, not across.

Sources: results/ldsc/h2.tsv (standard two-step LDSC) and the FinnGen R12 manifest. Every sentence in the
reading panel is computed.
"""
import numpy as np
import pandas as pd

import itoju_svg as sv
from itoju_labels import GROUPS

ROOT = sv.ROOT
SETS = ["Epilepsy", "Sleep apnoea", "Insomnia", "Constipation", "ADHD", "Intellectual disability"]
SHORT = {
    "G6_EPLEPSY": "any", "FE": "focal", "FE_STRICT": "focal, strict", "FE_MODE": "focal, mode",
    "GE": "generalized", "GE_STRICT": "generalized, strict", "GE_MODE": "generalized, mode",
    "G6_SLEEPAPNO": "hospital records", "G6_SLEEPAPNO_INCLAVO": "+ primary care", "SLEEP": "any sleep disorder",
    "F5_INSOMNIA": "F51.0 or G47.0", "KRA_PSY_SLEEP_NONORG_EXMORE": "all of F51",
    "K11_CONSTIPATION": "K59.0 or laxatives", "K11_OTHFUNC": "all of K59",
    "F5_ADHD": "F90.0", "KRA_PSY_HYPERKIN_EXMORE": "all of F90",
    "F5_MILDRET": "mild F70", "KRA_PSY_MENTALRET_EXMORE": "any (F7)",
}
YMAX = 0.25
Z50, Z95 = 0.674, 1.96
TICKS = [(1e3, "1k"), (2e3, "2k"), (5e3, "5k"), (1e4, "10k"), (2e4, "20k"), (5e4, "50k"), (1e5, "100k")]
COLW = [176.0, 158.0, 158.0]
PLOTW = [52.0, 40.0, 40.0]
ROW_TOP = [114.0, 268.0]
PH = 96.0
LABEL_GAP = 9.0


def main():
    h2 = pd.read_csv(ROOT / "results" / "ldsc" / "h2.tsv", sep="\t").set_index("trait")
    man = pd.read_csv(ROOT / "data" / "raw" / "finngen_R12_manifest.tsv", sep="\t").set_index("phenocode")
    groups = dict(GROUPS)
    n_cases = lambda e: int(man.loc[e, "num_cases"])
    val = lambda e: float(h2.loc[e, "h2"])

    W = sv.DOUBLE
    left = 14.0
    H = 490.0
    f = sv.Figure(W, H)
    f.header("itoju  /  heritability",
             "Heritability follows what a definition catches, not how many cases it has",
             "Liability-scale SNP heritability (up) against the number of cases (across, log scale), per definition")
    f.key_row(left, 70, [("dot", sv.FIRST, "first definition of the condition"),
                         ("dot", sv.SECOND, "alternative definitions"),
                         ("hollow", sv.DIM, "estimate below zero, drawn at zero")])
    f.key_row(left, 83, [("line", sv.INK_2, "bar: 50% interval"),
                         ("band", sv.NOISE, "sleeve: 95% interval, cut at zero")])

    for p, name in enumerate(SETS):
        c, row = p % 3, p // 3
        x0 = left + sum(COLW[:c])
        rt = ROW_TOP[row]
        defs = groups[name]
        cases = np.array([n_cases(e) for e, _ in defs], float)
        lc = np.log10(cases)
        mid = (lc.min() + lc.max()) / 2
        half = max(0.35, (lc.max() - lc.min()) / 2 + 0.12)
        px0, px1 = x0 + 22, x0 + 22 + PLOTW[c]
        py0, py1 = rt + 18, rt + 18 + PH
        X = sv.scale(mid - half, mid + half, px0, px1)
        Y = sv.scale(0, YMAX, py1, py0)

        f.text(x0, rt, name, 8.2, sv.INK, family=sv.SERIF)
        f.text(x0, rt + 10, f"{int(cases.min()):,} to {int(cases.max()):,} cases", 7.0, sv.DIM)
        for v in (0.1, 0.2):
            f.line(px0, Y(v), px1, Y(v), sv.GRID, 0.5)
        for v, s in ((0, "0"), (0.1, "0.1"), (0.2, "0.2")):
            f.text(px0 - 4, Y(v) + 2.4, s, 7.0, sv.DIM, anchor="end")
        f.line(px0, py1, px1, py1, sv.RULE, 0.6)
        last = -1e9
        for v, s in TICKS:
            xv = X(np.log10(v))
            if mid - half <= np.log10(v) <= mid + half and xv - last >= 16:     # skip a tick that would crowd the last
                f.line(xv, py1, xv, py1 + 2.5, sv.RULE, 0.5)
                f.text(xv, py1 + 10.5, s, 7.0, sv.DIM, anchor="middle")
                last = xv

        pts = []
        for i, (e, _) in enumerate(defs):
            r = h2.loc[e]
            col = sv.FIRST if i == 0 else sv.SECOND
            x = X(lc[i])
            lo95, hi95 = max(r.h2 - Z95 * r.h2_se, 0), min(r.h2 + Z95 * r.h2_se, YMAX)
            f.rect(x - 2.6, Y(hi95), 5.2, Y(lo95) - Y(hi95), sv.NOISE, rx=1, mark=True)
        for i, (e, _) in enumerate(defs):
            r = h2.loc[e]
            col = sv.FIRST if i == 0 else sv.SECOND
            x = X(lc[i])
            if r.h2 > 0:
                f.line(x, Y(max(r.h2 - Z50 * r.h2_se, 0)), x, Y(r.h2 + Z50 * r.h2_se), col, 2.0, mark=True)
                f.circle(x, Y(r.h2), 2.5, col, stroke=sv.GROUND, sw=0.6)
                pts.append((Y(r.h2), SHORT[e], col))
            else:
                f.circle(x, Y(0), 2.3, sv.GROUND, stroke=sv.DIM, sw=0.9)
                pts.append((Y(0), SHORT[e] + ", below 0", sv.DIM))

        # labels stacked beside the panel, nudged apart only as far as needed; the elbow starts at the dot's height
        pts.sort()
        ys = []
        for yv, _, _ in pts:
            ys.append(max(yv, ys[-1] + LABEL_GAP) if ys else yv)
        over = ys[-1] - (py1 + 1)
        if over > 0:
            ys = [v - over for v in ys]
        lx = px1 + 10
        for (yv, lab, col), yl in zip(pts, ys):
            f.polyline([(px1 + 2, yv), (px1 + 5, yv), (lx - 2, yl)], sv.RULE, 0.5, mark=False)
            f.circle(lx + 1.6, yl - 0.1, 1.6, col, mark=False)
            f.text(lx + 6, yl + 2.4, lab, 7.0, sv.INK_2, max_w=COLW[c] - (lx + 6 - x0) - 2)

    lower = []
    for name in SETS:
        es = [e for e, _ in groups[name]]
        most, fewest = max(es, key=n_cases), min(es, key=n_cases)
        if val(most) < val(fewest):
            lower.append(name.lower())
    print(f"  most cases, lower h2: {lower}; GE {val('GE'):.4f} FE {val('FE'):.4f} any {val('G6_EPLEPSY'):.4f}")
    lines = [
        f"Sleep apnoea: {val('G6_SLEEPAPNO'):.3f} from hospital records, {val('G6_SLEEPAPNO_INCLAVO'):.3f} with primary care, "
        f"{val('SLEEP'):.3f} for any sleep disorder.",
        f"Epilepsy: generalized {val('GE'):.3f} on {n_cases('GE'):,} cases, focal {val('FE'):.3f} on {n_cases('FE'):,}; "
        f"any epilepsy {val('G6_EPLEPSY'):.3f} sits between.",
        f"Most cases, lower heritability: {len(lower)} of {len(SETS)} conditions ({', '.join(lower)}).",
        "Sample prevalence stands in for population prevalence, so compare heights within a condition.",
    ]
    f.reading(left, 414, W - left - 10, lines, strong=(0,))
    f.source("results/ldsc/h2.tsv, two-step LDSC; case counts from the FinnGen R12 manifest")
    f.save("fig08_heritability")


if __name__ == "__main__":
    main()
