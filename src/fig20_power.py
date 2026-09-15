"""Figure: what buys the power to see a definition effect. Shared cases, not sample size.

The smallest difference between two genetic correlations detectable at alpha 0.05 with 80% power is
    MDD = (z_0.975 + z_0.80) * sqrt(se1^2 + se2^2 - 2 * rho * se1 * se2),
where rho is the correlation between the two estimates' sampling errors. The shaded ground is that formula
for two estimates with the same SE (vertical axis, log scale) and error correlation rho (horizontal), drawn
as a contour map: each line is one detectable difference. Definitions that share most of their cases sit
far right, where the contours plunge.

Each dot is one comparison between two definitions of the same condition (135), in the colour of its
psychiatric trait, at its jackknife error correlation and the root mean square of its two SEs; turmeric
rings are the 31 EHR against clinical depression comparisons. With unequal SEs the map is an approximation,
so the script reports how many dots sit in the band of their exact value. Four comparisons are numbered.

Sources: results/definition_effect/pairwise_delta.tsv and rg_common.tsv.
"""
import numpy as np
import pandas as pd
from scipy.stats import norm

import itoju_svg as sv

ROOT = sv.ROOT
K = norm.ppf(0.975) + norm.ppf(0.80)
LEVELS = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 10.0]
BANDS = ["#E9EFEC", "#DBE4E0", "#CBD7D2", "#B7C6C0", "#9FB2AB", "#839890", "#677D74"]
RHO = (-0.2, 1.0)
SE = (0.012, 0.45)
EXAMPLES = [  # trait, def1, def2, name
    ("MDD", "G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "sleep apnoea, adding primary care (depression)"),
    ("PTSD", "K11_CONSTIPATION", "K11_OTHFUNC", "constipation or laxatives vs all of K59 (PTSD)"),
    ("MDD", "FE", "GE", "focal vs generalized epilepsy (depression)"),
    ("PTSD", "MDD_EHR", "MDD_Clin", "EHR vs clinical depression (PTSD)"),
]


def main():
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    rc = pd.read_csv(ROOT / "results" / "definition_effect" / "rg_common.tsv", sep="\t")
    se = {}
    for r in rc.itertuples():
        se[(r.trait1, r.trait2)] = se[(r.trait2, r.trait1)] = r.se
    pw["se1"] = [se[(r.shared, r.def1)] for r in pw.itertuples()]
    pw["se2"] = [se[(r.shared, r.def2)] for r in pw.itertuples()]
    pw["se_rms"] = np.sqrt((pw.se1 ** 2 + pw.se2 ** 2) / 2)
    exact = K * np.sqrt(pw.se1 ** 2 + pw.se2 ** 2 - 2 * pw.error_corr * pw.se1 * pw.se2)
    assert np.max(np.abs(exact - pw.min_detectable_delta)) < 1e-3
    band_of = lambda v: np.searchsorted(LEVELS, v) - 1
    same_band = float(np.mean(band_of(exact.values) == band_of(K * pw.se_rms.values * np.sqrt(2 - 2 * pw.error_corr.values))))
    ctrl = pw.family == "MDD definition (positive control)"

    W = sv.DOUBLE
    left = 14.0
    px0, px1, py0, py1 = 48.0, 318.0, 74.0, 344.0
    X = sv.scale(RHO[0], RHO[1], px0, px1)
    ly = lambda v: np.log10(v)
    Yl = sv.scale(ly(SE[1]), ly(SE[0]), py0, py1)
    Y = lambda v: Yl(ly(v))
    H = 448.0
    f = sv.Figure(W, H)
    f.header("itoju  /  what buys power",
             "Shared cases, not sample size, decide how small a shift can be seen",
             "Smallest difference in genetic correlation detectable 4 times in 5, by error correlation and SE")

    # equal SEs: detectable = K * s * sqrt(2 - 2 rho), so each contour is s = L / (K sqrt(2 - 2 rho)), drawn exactly
    rr = np.concatenate([np.linspace(RHO[0], 0.9, 400), 1 - np.logspace(-1, -7, 400)[1:]])
    s_of = lambda lv: lv / (K * np.sqrt(2 - 2 * rr))
    clip = lambda v: np.clip(v, SE[0], SE[1])
    for i in range(len(LEVELS) - 1):
        lo = clip(s_of(LEVELS[i])) if LEVELS[i] > 0 else np.full_like(rr, SE[0])
        hi = clip(s_of(LEVELS[i + 1]))
        pts = [(X(a), Y(b)) for a, b in zip(rr, lo)] + [(X(a), Y(b)) for a, b in zip(rr[::-1], hi[::-1])]
        f.polygon(pts, BANDS[i])
    for lv in LEVELS[1:-1]:
        sv_ = s_of(lv)
        m = (sv_ >= SE[0]) & (sv_ <= SE[1])
        f.polyline([(X(a), Y(b)) for a, b in zip(rr[m], sv_[m])], "#FFFFFF", 0.8, mark=False, opacity=0.95)
    f.rect(px0, py0, px1 - px0, py1 - py0, "none", stroke=sv.RULE, sw=0.5)
    x_lab = RHO[0] + 0.012
    for lv in LEVELS[1:-1]:
        y_lab = lv / (K * np.sqrt(2 - 2 * x_lab))
        ink = sv.INK if lv <= 0.2 else "#FFFFFF"
        if SE[0] * 1.2 < y_lab < SE[1] * 0.85:
            f.text(X(x_lab) + 1, Y(y_lab) - 2, f"{lv:g}", 7.0, ink)
        else:
            r_lab = 1 - (lv / (K * SE[0] * 1.15)) ** 2 / 2
            if RHO[0] + 0.05 < r_lab < 0.95:
                f.text(X(r_lab) + 2, Y(SE[0] * 1.15) - 2, f"{lv:g}", 7.0, ink)

    for r in pw[~ctrl].itertuples():
        f.circle(X(r.error_corr), Y(r.se_rms), 1.9, sv.TRAIT[r.shared], stroke="#FFFFFF", sw=0.35)
    for r in pw[ctrl].itertuples():
        f.circle(X(r.error_corr), Y(r.se_rms), 1.7, "none", stroke=sv.TRAIT["MDD"], sw=0.9)
    ex_rows = []
    for k, (t, d1, d2, name) in enumerate(EXAMPLES, start=1):
        r = pw[(pw.shared == t) & (pw.def1 == d1) & (pw.def2 == d2)].iloc[0]
        cx, cy = X(r.error_corr), Y(r.se_rms)
        f.circle(cx, cy, 4.2, "none", stroke=sv.INK, sw=0.8, mark=False)
        dx = -9 if r.error_corr > 0.8 else 9
        f.circle(cx + dx, cy - 8, 4.3, "#FFFFFF", stroke=sv.INK, sw=0.6, mark=False)
        f.text(cx + dx, cy - 5.6, str(k), 7.0, sv.INK, anchor="middle", on_mark=True)
        ex_rows.append((k, name, r.error_corr, r.se_rms, r.min_detectable_delta))
        print(f"  example {k}: rho {r.error_corr:.3f}, se_rms {r.se_rms:.4f}, detectable {r.min_detectable_delta:.3f}")

    for v in (-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0):
        f.line(X(v), py1, X(v), py1 + 2.5, sv.RULE, 0.5)
        f.text(X(v), py1 + 10.5, f"{v:g}", 7.0, sv.DIM, anchor="middle")
    f.text((px0 + px1) / 2, py1 + 21, "error correlation between the two estimates (shared cases push it right)", 7.0,
           sv.DIM, anchor="middle")
    for v in (0.02, 0.05, 0.1, 0.2, 0.4):
        f.line(px0 - 2.5, Y(v), px0, Y(v), sv.RULE, 0.5)
        f.text(px0 - 4, Y(v) + 2.4, f"{v:g}", 7.0, sv.DIM, anchor="end")
    f.text(px0 - 4, py0 - 6, "SE of each estimate", 7.0, sv.DIM)

    rx = 336.0
    f.key_traits(rx, py0 + 4, ["ASD", "SCZ", "BIP"])
    f.key_traits(rx, py0 + 16, ["MDD", "PTSD"])
    f.key_row(rx, py0 + 28, [("hollow", sv.TRAIT["MDD"], "EHR vs clinical depression")])
    f.key_row(rx, py0 + 40, [("band", BANDS[2], "ground: detectable difference")])
    y = py0 + 64
    f.text(rx, y, "NUMBERED COMPARISONS", 7.0, sv.DIM, spacing=0.6)
    y += 13
    for k, name, rho, s, mdd in ex_rows:
        words, line_, out = name.split(), "", []
        for w_ in words:
            if sv.text_width((line_ + " " + w_).strip()) > W - 10 - rx - 14:
                out.append(line_)
                line_ = w_
            else:
                line_ = (line_ + " " + w_).strip()
        out.append(line_)
        f.circle(rx + 4, y - 2.5, 4.3, "#FFFFFF", stroke=sv.INK, sw=0.6, mark=False)
        f.text(rx + 4, y, str(k), 7.0, sv.INK, anchor="middle", on_mark=True)
        for i, s_ in enumerate(out):
            f.text(rx + 14, y + i * 9, s_, 7.0, sv.INK_2)
        f.text(rx + 14, y + len(out) * 9, f"detects {mdd:.3f}, error corr {rho:.2f}" if mdd < 0.1 else f"detects {mdd:.2f}, error corr {rho:.2f}", 7.0, sv.INK)
        y += (len(out) + 1) * 9 + 8

    lines = [
        f"The same SE buys a shift of {K * 0.05 * np.sqrt(2):.2f} with no shared cases, and {K * 0.05 * np.sqrt(2 - 2 * 0.95):.2f} "
        "when the errors correlate 0.95.",
        f"Sleep apnoea, adding primary care: SE {ex_rows[0][3]:.3f}, detects {ex_rows[0][4]:.3f}. "
        f"EHR vs clinical: SE {ex_rows[3][3]:.3f}, detects {ex_rows[3][4]:.2f}.",
        f"The map assumes equal SEs; {same_band:.0%} of comparisons sit in the band of their exact detectable shift.",
    ]
    f.reading(left, 378, W - left - 10, lines, strong=(0,))
    f.source("results/definition_effect/pairwise_delta.tsv; detectable = (z0.975 + z0.80) x SE of the difference")
    f.save("fig20_power")
    paper_single(pw, ctrl)


def paper_single(pw, ctrl):
    """The paper's single-column layout: the key above, the same map, the numbered comparisons below."""
    W = sv.SINGLE
    left, right = 10.0, sv.SINGLE - 8
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    f.key_traits(left, 12, ["ASD", "SCZ", "BIP"])
    kx = f.key_traits(left, 24, ["MDD", "PTSD"])
    f.text(kx, 24, "white lines: detectable shift", 7.0, sv.INK_2, max_w=right - kx)
    f.key_row(left, 36, [("hollow", sv.TRAIT["MDD"], "EHR vs clinical depression")])
    px0, px1, py0, py1 = 36.0, right, 58.0, 258.0
    X = sv.scale(RHO[0], RHO[1], px0, px1)
    Yl = sv.scale(np.log10(SE[1]), np.log10(SE[0]), py0, py1)
    Y = lambda v: Yl(np.log10(v))
    rr = np.concatenate([np.linspace(RHO[0], 0.9, 300), 1 - np.logspace(-1, -7, 300)[1:]])
    s_of = lambda lv: lv / (K * np.sqrt(2 - 2 * rr))
    clip = lambda v: np.clip(v, SE[0], SE[1])
    for i in range(len(LEVELS) - 1):
        lo = clip(s_of(LEVELS[i])) if LEVELS[i] > 0 else np.full_like(rr, SE[0])
        hi = clip(s_of(LEVELS[i + 1]))
        f.polygon([(X(a), Y(b)) for a, b in zip(rr, lo)] + [(X(a), Y(b)) for a, b in zip(rr[::-1], hi[::-1])], BANDS[i])
    for lv in LEVELS[1:-1]:
        sv_ = s_of(lv)
        m = (sv_ >= SE[0]) & (sv_ <= SE[1])
        f.polyline([(X(a), Y(b)) for a, b in zip(rr[m], sv_[m])], "#FFFFFF", 0.7, mark=False, opacity=0.95)
    f.rect(px0, py0, px1 - px0, py1 - py0, "none", stroke=sv.RULE, sw=0.5)
    x_lab = RHO[0] + 0.015
    for lv in LEVELS[1:-1]:
        y_lab = lv / (K * np.sqrt(2 - 2 * x_lab))
        ink = sv.INK if lv <= 0.2 else "#FFFFFF"
        if SE[0] * 1.25 < y_lab < SE[1] * 0.8:
            f.text(X(x_lab) + 1, Y(y_lab) - 2, f"{lv:g}", 7.0, ink)
    for r in pw[~ctrl].itertuples():
        f.circle(X(r.error_corr), Y(r.se_rms), 1.5, sv.TRAIT[r.shared], stroke="#FFFFFF", sw=0.3)
    for r in pw[ctrl].itertuples():
        f.circle(X(r.error_corr), Y(r.se_rms), 1.4, "none", stroke=sv.TRAIT["MDD"], sw=0.8)
    ex = []
    for k, (t, d1, d2, name) in enumerate(EXAMPLES, start=1):
        r = pw[(pw.shared == t) & (pw.def1 == d1) & (pw.def2 == d2)].iloc[0]
        cx, cy = X(r.error_corr), Y(r.se_rms)
        f.circle(cx, cy, 3.6, "none", stroke=sv.INK, sw=0.7, mark=False)
        dx = -8 if r.error_corr > 0.8 else 8
        f.circle(cx + dx, cy - 7, 4.1, "#FFFFFF", stroke=sv.INK, sw=0.6, mark=False)
        f.text(cx + dx, cy - 4.6, str(k), 7.0, sv.INK, anchor="middle", on_mark=True)
        ex.append((k, name, r.error_corr, r.min_detectable_delta))
    for v in (-0.2, 0.2, 0.6, 1.0):
        f.line(X(v), py1, X(v), py1 + 2.5, sv.RULE, 0.5)
        f.text(X(v), py1 + 10.5, f"{v:g}", 7.0, sv.DIM, anchor="middle")
    f.text((px0 + px1) / 2, py1 + 20, "error correlation (shared cases push it right)", 7.0, sv.DIM, anchor="middle")
    for v in (0.02, 0.05, 0.1, 0.2, 0.4):
        f.line(px0 - 2.5, Y(v), px0, Y(v), sv.RULE, 0.5)
        f.text(px0 - 4, Y(v) + 2.4, f"{v:g}", 7.0, sv.DIM, anchor="end")
    f.text(px0 - 4, py0 - 5, "SE of each estimate", 7.0, sv.DIM)
    y = py1 + 38
    for k, name, rho, mdd in ex:
        f.circle(left + 4, y - 2.5, 4.1, "#FFFFFF", stroke=sv.INK, sw=0.6, mark=False)
        f.text(left + 4, y, str(k), 7.0, sv.INK, anchor="middle", on_mark=True)
        f.text(left + 13, y, name, 7.0, sv.INK_2, max_w=right - left - 13)
        f.text(left + 13, y + 9, f"detects {mdd:.3f}, error corr {rho:.2f}" if mdd < 0.1 else f"detects {mdd:.2f}, error corr {rho:.2f}",
               7.0, sv.INK, max_w=right - left - 13)
        y += 22
    f.h = y - 10
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.save("fig20_power_paper", paper=False, png=False)


if __name__ == "__main__":
    main()
