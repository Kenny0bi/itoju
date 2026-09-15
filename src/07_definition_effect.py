"""The core test: does the genetic correlation move when only the case definition changes?

Design, fixed in docs/ANALYSIS_PLAN.md before any rg was estimated:
  DEFINITION_SETS  alternative definitions of one construct (they overlap or nest)
  NEIGHBOURS       different conditions, kept as a comparison for how far rg moves between conditions
  SHARED           the trait each definition is correlated with (autism is focal)
  POSITIVE CONTROL MDD_EHR vs MDD_Clin, two definitions of depression with a known definition effect

Steps
  1. Restrict every trait to one common SNP set (non-missing Z in all), so LDSC's 200 jackknife
     blocks are identical across pairs.
  2. For each shared trait X, run one LDSC call X,D1,D2,... with --two-step 99999 and
     --print-delete-vals (one-step estimator, evenly spaced blocks).
  3. Rebuild per-block rg pseudovalues exactly as LDSC's RatioJackknife does, and check that they
     reproduce LDSC's printed SE.
  4. Within each definition set and shared trait:
       - pairwise delta rg, SE from the pseudovalue difference, z, p, error correlation,
         minimum detectable delta (alpha 0.05, power 0.80);
       - heterogeneity Q = r' W r - (1' W r)^2 / (1' W 1), W = inverse jackknife covariance,
         df = k - 1 (pseudo-inverse if the covariance is near singular; condition number reported).
  5. Join pairwise deltas to the definition distances from src/02_definition_anatomy.py.

Outputs (results/definition_effect/): rg_common.tsv, pairwise_delta.tsv, heterogeneity.tsv
"""
import os
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT.parent / "tools"
PY = TOOLS / "ldsc_env" / "bin" / "python"
LDSC = TOOLS / "cbiit_ldsc" / "ldsc.py"
LD = str(ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr") + "/"
MUNGED = ROOT / "data" / "munged"
COMMON = ROOT / "data" / "munged_common"
OUT = ROOT / "results" / "definition_effect"
RUNS = OUT / "ldsc"

DEFINITION_SETS = {
    "Epilepsy": ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"],
    "Sleep apnoea": ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"],
    "Insomnia": ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"],
    "Constipation": ["K11_CONSTIPATION", "K11_OTHFUNC"],
    "ADHD": ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"],
    "Intellectual disability": ["F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE"],
}
NEIGHBOURS = ["G6_STATUSEPI", "G6_SLEEPDISOTH", "F5_SLEEP_NOS", "K11_IBS", "K11_FUNCDYSP", "K11_REFLUX",
              "N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD", "KRA_PSY_DEVWIDE_EXMORE", "KRA_PSY_AUTISM_EXMORE"]
SHARED = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
CONTROL = ["MDD_EHR", "MDD_Clin"]
FILES = {"ASD": "ASD_Grove2019"}

ENV = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2",
           VECLIB_MAXIMUM_THREADS="2")
Z_MDD = norm.ppf(0.975) + norm.ppf(0.80)


def src(t):
    return MUNGED / f"{FILES.get(t, t)}.sumstats.gz"


def build_common(traits):
    """Write copies of every trait restricted to SNPs with non-missing Z in all of them."""
    COMMON.mkdir(parents=True, exist_ok=True)
    marker = COMMON / "snps.txt"
    frames = {t: pd.read_csv(src(t), sep="\t").dropna(subset=["Z"]) for t in traits}
    common = set.intersection(*(set(f.SNP) for f in frames.values()))
    if marker.exists() and set(marker.read_text().split()) == common:
        return len(common)
    for t, f in frames.items():
        f[f.SNP.isin(common)].to_csv(COMMON / f"{t}.sumstats.gz", sep="\t", index=False, compression="gzip")
    marker.write_text("\n".join(sorted(common)))
    return len(common)


def run_rg(first, others, tag=""):
    """One LDSC call: first against each of others, one-step, delete values printed.
    `tag` keeps runs with the same first trait but different partners from sharing an output name
    (the positive-control runs would otherwise collide with the endpoint runs)."""
    RUNS.mkdir(parents=True, exist_ok=True)
    out = RUNS / f"rg_{first}{tag}"
    log = Path(str(out) + ".log")
    if log.exists() and log.read_text().count("Genetic Correlation: ") >= len(others) \
            and "Summary of Genetic Correlation Results" in log.read_text():
        return log
    paths = ",".join(str(COMMON / f"{t}.sumstats.gz") for t in [first, *others])
    cmd = ["nice", "-n", "10", str(PY), str(LDSC), "--rg", paths, "--ref-ld-chr", LD, "--w-ld-chr", LD,
           "--two-step", "99999", "--print-delete-vals", "--out", str(out)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ENV)
    return log


def parse_run(first, others, log):
    """Summary table plus pseudovalues for each pair in one LDSC run."""
    text = log.read_text()
    lines = text[text.find("Summary of Genetic Correlation Results"):].splitlines()[1:]
    hdr = lines[0].split()
    table = [dict(zip(hdr, l.split())) for l in lines[1:1 + len(others)]]
    rows, pseudo = [], {}
    for t, rec in zip(others, table):
        stem = f"{log.stem}{first}.sumstats.gz_{t}.sumstats.gz"
        rg = float(rec["rg"]) if rec["rg"] != "NA" else np.nan
        se = float(rec["se"]) if rec["se"] != "NA" else np.nan
        row = {"trait1": first, "trait2": t, "rg": rg, "se": se, "p": float(rec["p"]) if rec["p"] != "NA" else np.nan,
               "h2_obs_2": float(rec["h2_obs"]), "gcov_int": float(rec["gcov_int"]),
               "gcov_int_se": float(rec["gcov_int_se"])}
        if np.isfinite(rg):
            g = np.loadtxt(RUNS / f"{stem}.gencov.delete").reshape(-1)
            h1 = np.loadtxt(RUNS / f"{stem}.hsq1.delete").reshape(-1)
            h2 = np.loadtxt(RUNS / f"{stem}.hsq2.delete").reshape(-1)
            n = len(g)
            with np.errstate(invalid="ignore"):
                pv = n * rg - (n - 1) * g / np.sqrt(h1 * h2)
            rebuilt = float(np.sqrt(np.var(pv, ddof=1) / n))
            row["se_rebuilt"] = rebuilt
            row["se_check_ok"] = bool(np.isfinite(rebuilt) and (abs(rebuilt - se) <= 6e-5 or abs(rebuilt - se) / se < 1e-3))
            pseudo[t] = pv
        rows.append(row)
    return rows, pseudo


def pairwise(x, family, defs, rg, pseudo):
    out = []
    for a, b in combinations(defs, 2):
        if a not in pseudo or b not in pseudo:
            continue
        d = rg[a] - rg[b]
        diff = pseudo[a] - pseudo[b]
        n = len(diff)
        se = float(np.sqrt(np.var(diff, ddof=1) / n))
        out.append({"shared": x, "family": family, "def1": a, "def2": b, "rg1": rg[a], "rg2": rg[b],
                    "delta": d, "se_delta": se, "z": d / se, "p": float(2 * norm.sf(abs(d / se))),
                    "error_corr": float(np.corrcoef(pseudo[a], pseudo[b])[0, 1]),
                    "min_detectable_delta": Z_MDD * se})
    return out


def heterogeneity(x, family, defs, rg, pseudo):
    defs = [d for d in defs if d in pseudo]
    if len(defs) < 2:
        return None
    r = np.array([rg[d] for d in defs])
    P = np.column_stack([pseudo[d] for d in defs])
    S = np.cov(P, rowvar=False, ddof=1) / P.shape[0]
    cond = float(np.linalg.cond(S))
    W = np.linalg.pinv(S)
    one = np.ones(len(r))
    pooled = float(one @ W @ r / (one @ W @ one))
    Q = float(r @ W @ r - (one @ W @ r) ** 2 / (one @ W @ one))
    df = len(r) - 1
    return {"shared": x, "family": family, "k": len(r), "pooled_rg": pooled, "Q": Q, "df": df,
            "p": float(chi2.sf(Q, df)), "rg_range": float(r.max() - r.min()), "cov_condition_number": cond}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    endpoints = [e for fam in DEFINITION_SETS.values() for e in fam] + NEIGHBOURS
    traits = [t for t in SHARED + CONTROL if src(t).exists()] + [e for e in endpoints if src(e).exists()]
    missing = [t for t in SHARED + CONTROL + endpoints if not src(t).exists()]
    if missing:
        print("not yet munged, skipped:", " ".join(missing))
    n_common = build_common(traits)
    print(f"common SNP set: {n_common:,} SNPs across {len(traits)} traits")

    present_ep = [e for e in endpoints if src(e).exists()]
    all_rows, pairs, het = [], [], []
    # shared psychiatric traits (and the two MDD definitions) against every FinnGen endpoint
    for x in [t for t in SHARED + CONTROL if src(t).exists()]:
        log = run_rg(x, present_ep)
        rows, pseudo = parse_run(x, present_ep, log)
        all_rows += rows
        rg = {r["trait2"]: r["rg"] for r in rows}
        if x in SHARED:
            for fam, defs in DEFINITION_SETS.items():
                pairs += pairwise(x, fam, defs, rg, pseudo)
                h = heterogeneity(x, fam, defs, rg, pseudo)
                if h:
                    het.append(h)

    # positive control: MDD_EHR vs MDD_Clin against each other trait (not MDD, which contains both)
    if all(src(t).exists() for t in CONTROL):
        others = [t for t in ["ASD", "SCZ", "BIP", "PTSD"] if src(t).exists()] + present_ep
        # FinnGen endpoints: reuse the MDD_EHR and MDD_Clin runs above (same common SNP set, one-step,
        # so identical blocks). Psychiatric traits: one call y,MDD_EHR,MDD_Clin.
        for y in others:
            rg_y, pv_y = {}, {}
            if y in present_ep:
                for d in CONTROL:
                    rows, pseudo = parse_run(d, present_ep, RUNS / f"rg_{d}.log")
                    rec = next(r for r in rows if r["trait2"] == y)
                    rg_y[d], pv_y[d] = rec["rg"], pseudo.get(y)
            else:
                log = run_rg(y, CONTROL, tag="__control")
                rows, pseudo = parse_run(y, CONTROL, log)
                all_rows += rows
                rg_y = {r["trait2"]: r["rg"] for r in rows}
                pv_y = pseudo
            if all(pv_y.get(d) is not None for d in CONTROL):
                pairs += pairwise(y, "MDD definition (positive control)", CONTROL, rg_y, pv_y)
        log = run_rg("MDD_EHR", ["MDD_Clin"], tag="__control")
        rows, _ = parse_run("MDD_EHR", ["MDD_Clin"], log)
        all_rows += rows

    rgc = pd.DataFrame(all_rows).drop_duplicates(["trait1", "trait2"])
    pw = pd.DataFrame(pairs)
    hz = pd.DataFrame(het)
    dist = pd.read_csv(ROOT / "results" / "definitions" / "pairwise_distance.tsv", sep="\t")
    rule = dist[["mode_differs", "conditions_differ", "controls_differ"]].any(axis=1).astype(float)
    dist["distance"] = 1 - (dist["icd10_jaccard"] + dist["source_jaccard"] + (1 - rule)) / 3
    fwd = dist.rename(columns={"def1": "def1", "def2": "def2"})
    rev = dist.rename(columns={"def1": "def2", "def2": "def1"})
    keep = ["def1", "def2", "distance", "icd10_jaccard", "source_jaccard", "log_case_ratio"]
    if not pw.empty:
        pw = pw.merge(pd.concat([fwd[keep], rev[keep]]).drop_duplicates(["def1", "def2"]),
                      on=["def1", "def2"], how="left")
    rgc.to_csv(OUT / "rg_common.tsv", sep="\t", index=False)
    pw.to_csv(OUT / "pairwise_delta.tsv", sep="\t", index=False)
    hz.to_csv(OUT / "heterogeneity.tsv", sep="\t", index=False)
    bad = rgc[rgc.get("se_check_ok", True) == False]  # noqa: E712
    if len(bad):
        raise SystemExit(f"{len(bad)} pairs failed the SE rebuild check:\n{bad}")
    pd.set_option("display.width", 250)
    print(hz.round(4).to_string())
    print(pw.round(4).to_string())


if __name__ == "__main__":
    main()
