"""Figure: why a difference between two genetic correlations can be precise.

LDSC estimates its standard errors with a jackknife over 200 genome blocks. On the common SNP set, block
k is the same stretch of genome for both definitions, so each block gives a pair of pseudovalues.
Top of each panel: the 200 pairs, centred on the two estimates, against the line of equal values.
Bottom: every block's difference (second definition minus first) as a tick, with its density above. The
coloured density keeps the blocks paired as they are; the gray one pairs each block with a random other block, which is what the
difference would look like if the two estimates shared no people. The narrower the filled swarm, the
smaller the standard error of the difference.

  a  depression: sleep apnoea (hospital records) against any sleep disorder, nested definitions
  b  PTSD: constipation or laxatives against all of K59, partly shared cases
  c  PTSD: depression from health records against clinical depression, almost no shared samples

Source: results/definition_effect/ldsc (LDSC delete values). The SE printed is the paper's estimator:
paired sqrt(var(p2 - p1) / n); pairing broken sqrt((var p1 + var p2) / n).
"""
import numpy as np

import itoju_svg as sv

ROOT = sv.ROOT
RUNS = ROOT / "results" / "definition_effect" / "ldsc"
PANELS = [
    ("a", "Nested definitions", "depression: sleep apnoea vs any sleep", "MDD", "G6_SLEEPAPNO", "SLEEP", "",
     "hospital records", "any sleep disorder"),
    ("b", "Partly shared cases", "PTSD: constipation vs all of K59", "PTSD", "K11_CONSTIPATION", "K11_OTHFUNC", "",
     "K59.0 or laxatives", "all of K59"),
    ("c", "Almost nothing shared", "PTSD: EHR vs clinical depression", "PTSD", "MDD_EHR", "MDD_Clin", "__control",
     "depression, EHR", "depression, clinical"),
]
LIM, DLIM, R_DOT = 2.4, 2.4, 1.05


def pseudo(run, second, tag=""):
    log = RUNS / f"rg_{run}{tag}.log"
    text = log.read_text()
    lines = text[text.find("Summary of Genetic Correlation Results"):].splitlines()[1:]
    hdr = lines[0].split()
    rec = next(dict(zip(hdr, l.split())) for l in lines[1:] if l.split() and l.split()[1].endswith(f"/{second}.sumstats.gz"))
    stem = f"{log.stem}{run}.sumstats.gz_{second}.sumstats.gz"
    g = np.loadtxt(RUNS / f"{stem}.gencov.delete").reshape(-1)
    h1 = np.loadtxt(RUNS / f"{stem}.hsq1.delete").reshape(-1)
    h2 = np.loadtxt(RUNS / f"{stem}.hsq2.delete").reshape(-1)
    n = len(g)
    rg = float(rec["rg"])
    return n * rg - (n - 1) * g / np.sqrt(h1 * h2), rg


def kde(values, grid):
    """Gaussian kernel density with Silverman's bandwidth."""
    bw = 1.06 * np.std(values, ddof=1) * len(values) ** (-0.2)
    z = (grid[:, None] - values[None, :]) / bw
    return np.exp(-0.5 * z * z).sum(axis=1) / (len(values) * bw * np.sqrt(2 * np.pi))


def swarm(xs, y0, d):
    placed, out = [], [None] * len(xs)
    for i in np.argsort(np.abs(xs - np.median(xs))):
        k = 0
        while True:
            off = ((k + 1) // 2) * d * (1 if k % 2 else -1)
            if all((xs[i] - px) ** 2 + (off - py) ** 2 >= d * d for px, py in placed):
                break
            k += 1
        placed.append((xs[i], off))
        out[i] = (xs[i], y0 + off)
    return out


def main():
    W = sv.DOUBLE
    left, gap = 14.0, 14.0
    col_w = (W - left - 10 - 2 * gap) / 3
    H = 522.0
    f = sv.Figure(W, H)
    f.header("itoju  /  why the difference is precise",
             "When two definitions share people, their noise cancels block by block",
             "LDSC's 200-block jackknife, paired block by block across two definitions of one condition")
    f.key_row(left, 72, [("dot", sv.INK_2, "one genome block, both definitions"),
                         ("line", sv.INK_2, "one block's difference"),
                         ("band", sv.NOISE, "blocks paired at random")])
    readings = []
    for k, (letter, title, sub, trait, d1, d2, tag, lab1, lab2) in enumerate(PANELS):
        x0 = left + k * (col_w + gap)
        col = sv.TRAIT[trait]
        p1, r1 = pseudo(trait, d1, tag)
        p2, r2 = pseudo(trait, d2, tag)
        n = len(p1)
        dev1, dev2 = p1 - r1, p2 - r2
        corr = float(np.corrcoef(dev1, dev2)[0, 1])
        se_pair = float(np.sqrt(np.var(p2 - p1, ddof=1) / n))
        se_broken = float(np.sqrt((np.var(p1, ddof=1) + np.var(p2, ddof=1)) / n))
        readings.append((title, sub, se_pair, se_broken, corr))
        print(f"  {letter} {trait} {d1} vs {d2}: r {corr:.3f}, SE paired {se_pair:.4f}, broken {se_broken:.4f}")

        f.text(x0, 96, f"{letter}  {title}", 9.0, sv.INK, family=sv.SERIF)
        f.text(x0, 106.5, sub, 7.0, sv.INK_2)
        side = 116.0
        sx0, sy0 = x0 + (col_w - side) / 2, 116.0
        S = sv.scale(-LIM, LIM, 0, side)
        f.rect(sx0, sy0, side, side, "none", stroke=sv.GRID, sw=0.6)
        f.line(sx0, sy0 + side, sx0 + side, sy0, sv.RULE, 0.6, dash="2 2")
        for a, b in zip(np.clip(dev1, -LIM, LIM), np.clip(dev2, -LIM, LIM)):
            f.circle(sx0 + S(a), sy0 + side - S(b), 1.15, col, opacity=0.85)
        f.text(sx0, sy0 + side + 10, f"across: {lab1}", 7.0, sv.DIM)
        f.text(sx0, sy0 + side + 19, f"up: {lab2}", 7.0, sv.DIM)
        off = int(((np.abs(dev1) > LIM) | (np.abs(dev2) > LIM)).sum())
        f.text(x0, sy0 + side + 33, f"blocks agree, r = {corr:.2f}", 7.2, sv.INK)
        if off:
            f.text(x0 + col_w, sy0 + side + 33, f"{off} off scale", 7.0, sv.DIM, anchor="end")

        D = sv.scale(-DLIM, DLIM, x0 + 4, x0 + col_w - 4)
        rng = np.random.default_rng(7)
        diff = dev2 - dev1
        broken = dev2[rng.permutation(n)] - dev1
        grid = np.linspace(-DLIM, DLIM, 241)
        dens = [kde(diff, grid), kde(broken, grid)]
        peak = max(d.max() for d in dens)
        y_pair, y_broken, y_axis = 334.0, 404.0, 420.0
        f.text(x0, 284, "paired", 7.0, sv.INK_2)
        f.text(x0 + col_w, 284, f"SE {se_pair:.3f}", 9.0, sv.INK, family=sv.SERIF, anchor="end")
        f.text(x0, 354, "pairing broken", 7.0, sv.INK_2)
        f.text(x0 + col_w, 354, f"SE {se_broken:.3f}", 9.0, sv.INK, family=sv.SERIF, anchor="end")
        for d, yb, fill, stroke in ((dens[0], y_pair, col, col), (dens[1], y_broken, sv.NOISE, sv.DIM)):
            h = 38.0 * d / peak
            pts = [(D(g), yb - hh) for g, hh in zip(grid, h)]
            f.polygon([(D(-DLIM), yb)] + pts + [(D(DLIM), yb)], fill, opacity=0.35 if fill == col else 1.0)
            f.polyline(pts, stroke, 0.9, mark=False)
        for v in np.clip(diff, -DLIM, DLIM):
            f.line(D(v), y_pair + 1.5, D(v), y_pair + 6, col, 0.35, opacity=0.8)
        for v in np.clip(broken, -DLIM, DLIM):
            f.line(D(v), y_broken + 1.5, D(v), y_broken + 6, sv.DIM, 0.35, opacity=0.8)
        f.line(x0 + 4, y_axis, x0 + col_w - 4, y_axis, sv.RULE, 0.5)
        for v, s in ((-2, "-2"), (0, "0"), (2, "+2")):
            f.line(D(v), y_axis, D(v), y_axis + 2.5, sv.RULE, 0.5)
            f.text(D(v), y_axis + 10.5, s, 7.0, sv.DIM, anchor="middle")
        f.line(D(0), 292, D(0), y_axis, sv.GRID, 0.5)
        f.text(x0 + col_w / 2, y_axis + 20, "per-block difference", 7.0, sv.DIM, anchor="middle")

    (t_a, _, sp_a, sb_a, _), (t_b, _, sp_b, sb_b, _), (t_c, _, sp_c, sb_c, _) = readings
    lines = [
        f"Nested definitions: pairing the blocks cuts the standard error of the difference from {sb_a:.3f} to {sp_a:.3f}.",
        f"Partly shared cases: from {sb_b:.3f} to {sp_b:.3f}. Almost nothing shared: from {sb_c:.3f} to {sp_c:.3f}.",
        "Shared people, not sample size, decide how small a definition effect can be seen.",
    ]
    f.reading(left, 452, W - left - 10, lines, strong=(0,))
    f.source("LDSC delete-one-block values on the common SNP set, 200 blocks")
    f.save("fig11_block_geometry")


if __name__ == "__main__":
    main()
