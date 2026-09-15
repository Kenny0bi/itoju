"""Real numbers for the jackknife part of the animation.

LDSC's standard error comes from a jackknife over 200 genome blocks. On the common SNP set, block k is
the same stretch of genome for every trait, so for two definitions of one condition the per-block
pseudovalues can be paired. This exports, for depression with sleep apnoea (hospital records) and with
any sleep disorder, the 200 pseudovalues under each definition, both estimates, the paired standard
error of the difference, and the standard error if the pairing is broken (as if the two estimates
shared no people). Values come from results/definition_effect/ldsc, the same files the paper uses.
Writes animation/block_numbers.json.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "results" / "definition_effect" / "ldsc"


def pseudo(first, second, tag=""):
    log = RUNS / f"rg_{first}{tag}.log"
    text = log.read_text()
    lines = text[text.find("Summary of Genetic Correlation Results"):].splitlines()[1:]
    hdr = lines[0].split()
    rec = None
    for l in lines[1:]:
        r = dict(zip(hdr, l.split()))
        if r.get("p2", "").endswith(f"/{second}.sumstats.gz"):
            rec = r
            break
    stem = f"{log.stem}{first}.sumstats.gz_{second}.sumstats.gz"
    g = np.loadtxt(RUNS / f"{stem}.gencov.delete").reshape(-1)
    h1 = np.loadtxt(RUNS / f"{stem}.hsq1.delete").reshape(-1)
    h2 = np.loadtxt(RUNS / f"{stem}.hsq2.delete").reshape(-1)
    n = len(g)
    rg = float(rec["rg"])
    return n * rg - (n - 1) * g / np.sqrt(h1 * h2), rg, float(rec["se"])


def main():
    p1, r1, s1 = pseudo("MDD", "G6_SLEEPAPNO")
    p2, r2, s2 = pseudo("MDD", "SLEEP")
    n = len(p1)
    se_pair = float(np.sqrt(np.var(p2 - p1, ddof=1) / n))
    se_broken = float(np.sqrt((np.var(p1, ddof=1) + np.var(p2, ddof=1)) / n))
    pw = pd.read_csv(ROOT / "results" / "definition_effect" / "pairwise_delta.tsv", sep="\t")
    row = pw[(pw.shared == "MDD") & (pw.def1 == "G6_SLEEPAPNO") & (pw.def2 == "SLEEP")].iloc[0]
    assert abs(se_pair - row.se_delta) < 1e-6, (se_pair, row.se_delta)
    assert abs((r2 - r1) + row.delta) < 1e-3
    out = dict(trait="MDD", def1="G6_SLEEPAPNO", def2="SLEEP", n_blocks=n,
               rg1=r1, rg2=r2, se1=s1, se2=s2, shift=r2 - r1, se_paired=se_pair, se_broken=se_broken,
               p=float(row.p), error_corr=float(row.error_corr), min_detectable=float(row.min_detectable_delta),
               dev1=(p1 - r1).round(4).tolist(), dev2=(p2 - r2).round(4).tolist())
    # which chromosome each jackknife block starts on: blocks are contiguous runs of the common SNP set in genome order
    snps = pd.read_csv(ROOT / "data" / "munged_common" / "MDD.sumstats.gz", sep="\t", usecols=["SNP"])
    ld = pd.concat([pd.read_csv(ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr" / f"{c}.l2.ldscore.gz", sep="\t",
                                usecols=["CHR", "SNP", "BP"]) for c in range(1, 23)])
    order = snps.merge(ld, on="SNP").sort_values(["CHR", "BP"]).reset_index(drop=True)
    m = len(order)
    starts = (np.arange(n) * m) // n
    out["n_snps_common"] = int(m)
    out["block_chrom"] = order.CHR.values[starts].astype(int).tolist()
    json.dump(out, open(ROOT / "animation" / "block_numbers.json", "w"), indent=1)
    print(f"blocks {n}; rg {r1:.4f} -> {r2:.4f}, shift {r2 - r1:+.4f}; SE paired {se_pair:.4f}, broken {se_broken:.4f}; "
          f"error corr {row.error_corr:.3f}; p {row.p:.2e}; pseudovalue sd {np.std(p1):.3f}, {np.std(p2):.3f}")


if __name__ == "__main__":
    main()
