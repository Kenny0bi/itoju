"""Standard error of a difference between two genetic correlations that share a trait.

Question: does rg(X, definition 1) differ from rg(X, definition 2) by more than chance allows?
The two estimates share trait X, and the two definitions share cases and controls, so their errors
are correlated. Treating them as independent would overstate the SE of the difference.

Method. LDSC estimates rg as gencov / sqrt(h2_1 * h2_2) and gets its SE from a 200-block
jackknife. When both pairs are run on one common SNP set with the one-step estimator
(--two-step 99999), their jackknife blocks are identical, so the block-level estimates line up:
  rg_del[k]  = gencov_del[k] / sqrt(h2_1_del[k] * h2_2_del[k])      (block k left out)
  pseudo[k]  = n * rg_hat - (n - 1) * rg_del[k]                    (LDSC RatioJackknife)
  delta[k]   = pseudo_A[k] - pseudo_B[k]
  SE(delta)  = sqrt(var(delta) / n)
The minimum detectable difference at two-sided alpha 0.05 and 80% power is
(z_0.975 + z_0.80) * SE(delta).

Self-check: before any difference is reported, the rebuilt single-pair SE must reproduce the SE
LDSC printed for that pair (relative gap below 1e-3), and the block counts must match.

Usage: python src/05_delta_rg.py LOG_A LOG_B
  where LOG_* are LDSC --rg logs run with --print-delete-vals on the common SNP set.
"""
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

Z_MDD = norm.ppf(0.975) + norm.ppf(0.80)


def load_pair(log_path):
    log_path = Path(log_path)
    text = log_path.read_text()
    tail = text[text.find("Summary of Genetic Correlation Results"):].splitlines()
    rec = dict(zip(tail[1].split(), tail[2].split()))
    stem = str(log_path)[:-len(".log")]
    deletes = {}
    for part in ("hsq1", "hsq2", "gencov"):
        hits = list(log_path.parent.glob(Path(stem).name + "*." + part + ".delete"))
        if len(hits) != 1:
            raise FileNotFoundError(f"expected one {part}.delete for {log_path.name}, found {len(hits)}")
        deletes[part] = np.loadtxt(hits[0]).reshape(-1)
    return {"rg": float(rec["rg"]), "se": float(rec["se"]), "p2": rec["p2"], **deletes}


def pseudovalues(pair):
    n = len(pair["gencov"])
    rg_del = pair["gencov"] / np.sqrt(pair["hsq1"] * pair["hsq2"])
    return n * pair["rg"] - (n - 1) * rg_del


def jackknife_se(pseudo):
    n = len(pseudo)
    return float(np.sqrt(np.var(pseudo, ddof=1) / n))


def compare(log_a, log_b):
    a, b = load_pair(log_a), load_pair(log_b)
    if len(a["gencov"]) != len(b["gencov"]):
        raise ValueError("block counts differ; the pairs were not run on the same SNP set")
    pa, pb = pseudovalues(a), pseudovalues(b)
    for name, pair, pseudo in (("A", a, pa), ("B", b, pb)):
        rebuilt = jackknife_se(pseudo)
        gap = abs(rebuilt - pair["se"]) / pair["se"]
        # LDSC prints SE to 4 decimals, so allow for rounding in the printed value
        if abs(rebuilt - pair["se"]) > 6e-5 and gap > 1e-3:
            raise ValueError(f"pair {name}: rebuilt SE {rebuilt:.5f} does not match LDSC SE {pair['se']:.4f}")
    delta = a["rg"] - b["rg"]
    se_delta = jackknife_se(pa - pb)
    corr = float(np.corrcoef(pa, pb)[0, 1])
    se_indep = float(np.sqrt(a["se"] ** 2 + b["se"] ** 2))
    return {
        "rg_a": a["rg"], "se_a": a["se"], "rg_b": b["rg"], "se_b": b["se"],
        "delta": delta, "se_delta": se_delta, "z": delta / se_delta,
        "p": float(2 * norm.sf(abs(delta / se_delta))),
        "error_corr": corr, "se_delta_if_independent": se_indep,
        "min_detectable_delta": Z_MDD * se_delta,
        "n_blocks": len(pa), "rebuilt_se_a": jackknife_se(pa), "rebuilt_se_b": jackknife_se(pb),
    }


if __name__ == "__main__":
    out = compare(sys.argv[1], sys.argv[2])
    for k, v in out.items():
        print(f"{k:26s} {v:.5f}" if isinstance(v, float) else f"{k:26s} {v}")
