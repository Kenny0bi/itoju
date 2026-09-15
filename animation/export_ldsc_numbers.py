"""Real numbers for the LD score regression animation.

The animation shows the one idea itoju rests on: for two GWAS, the product of a SNP's z-scores grows
with its LD score, the slope of that line is the genetic covariance, and the intercept is what
shared samples add. Two pairs with the same FinnGen endpoint (constipation) make the contrast:
  autism (iPSYCH-PGC, Danish)            x constipation: no shared samples, intercept near 0
  depression, EHR (includes FinnGen R5)  x constipation: shared samples, intercept above 0
SNPs on the common set are sorted into 20 equal-count LD score bins; each bin's mean z1*z2 is a
point. The straight lines are LDSC's own estimates from results/ldsc (two-step, full SNP sets):
  E[z1 z2 | l] = sqrt(N1 N2) * rho_g / M * l + intercept
so the animation draws the method's answer, not a refit of the bins.
Writes animation/ldsc_numbers.json.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "data" / "munged_common"
LD = ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr"


def ldscores():
    frames = [pd.read_csv(LD / f"{c}.l2.ldscore.gz", sep="\t", usecols=["SNP", "L2"]) for c in range(1, 23)]
    m = sum(float(open(LD / f"{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
    return pd.concat(frames), m


def rg_row(t1, t2):
    rg = pd.read_csv(ROOT / "results" / "ldsc" / "rg.tsv", sep="\t")
    return rg[(rg.trait1 == t1) & (rg.trait2 == t2)].iloc[0]


def gencov_obs(t1, t2):
    t = (ROOT / "results" / "ldsc" / f"rg_{t1}_{t2}.log").read_text()
    m = re.search(r"Total Observed scale gencov:\s*(-?[\d.]+)\s*\(([\d.]+)\)", t)
    return float(m.group(1)), float(m.group(2))


def main():
    l2, M = ldscores()
    endpoint = "K11_CONSTIPATION"
    out = {"endpoint": endpoint, "M_5_50": M, "pairs": []}
    for trait in ("ASD", "MDD_EHR"):
        a = pd.read_csv(COMMON / f"{trait}.sumstats.gz", sep="\t")
        b = pd.read_csv(COMMON / f"{endpoint}.sumstats.gz", sep="\t")
        m = a.merge(b, on="SNP", suffixes=("1", "2")).merge(l2, on="SNP")
        # align alleles: flip z2 where A1 differs
        flip = m.A11 != m.A12
        m.loc[flip, "Z2"] = -m.loc[flip, "Z2"]
        m["zz"] = m.Z1 * m.Z2
        m["bin"] = pd.qcut(m.L2.rank(method="first"), 20, labels=False)
        bins = m.groupby("bin").agg(l2=("L2", "mean"), zz=("zz", "mean"), n=("zz", "size")).reset_index()
        sem = m.groupby("bin").zz.std() / np.sqrt(bins.n)
        r = rg_row(trait, endpoint)
        g, gse = gencov_obs(trait, endpoint)
        n1, n2 = float(m.N1.median()), float(m.N2.median())
        slope = np.sqrt(n1 * n2) * g / M
        out["pairs"].append({
            "trait": trait, "n_snps": int(len(m)), "N1": n1, "N2": n2,
            "bins_l2": bins.l2.round(3).tolist(), "bins_zz": bins.zz.round(5).tolist(),
            "bins_sem": sem.round(5).tolist(),
            "gencov_obs": g, "gencov_obs_se": gse, "slope": slope, "intercept": float(r.gcov_int),
            "intercept_se": float(r.gcov_int_se), "rg": float(r.rg), "rg_se": float(r.se),
            "h2_obs_2": float(r.h2_obs),
        })
        print(f"{trait} x {endpoint}: {len(m):,} SNPs, slope {slope:.2e} per LD unit, intercept {r.gcov_int:.4f} "
              f"({r.gcov_int_se:.4f}), gencov {g:.4f} ({gse:.4f}), rg {r.rg:.3f} ({r.se:.3f}); "
              f"bin mean z1z2 {bins.zz.min():.4f} to {bins.zz.max():.4f}; LD {bins.l2.min():.1f} to {bins.l2.max():.1f}")
    json.dump(out, open(ROOT / "animation" / "ldsc_numbers.json", "w"), indent=1)


if __name__ == "__main__":
    main()
