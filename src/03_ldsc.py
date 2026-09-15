"""SNP heritability of every munged trait, and genetic correlations for the study design.

Traits:
  - psychiatric GWAS (focal ASD, plus SCZ, BIP, MDD, PTSD) and the MDD EHR and clinical subsets;
  - every munged FinnGen R12 endpoint.
Correlations:
  - each psychiatric trait x each FinnGen endpoint;
  - every pair of psychiatric traits (includes the positive control MDD_EHR x MDD_Clin).

Resumable: a run is skipped when its log already contains the result line. "Analysis finished" is
not enough, because the CBIIT port writes it even after a traceback.

Scales. rg runs carry no prevalence flags: rg is the same on either scale, and the CBIIT port
crashes in _get_rg_table when --pop-prev is given with --rg (NameError, `i` for `it`). Liability h2
comes from separate --h2 runs. FinnGen endpoints use K = sample prevalence, a stated approximation
for a biobank drawn largely from the population and hospital patients. PGC traits use effective N,
so their sample prevalence is 0.5, with population prevalences from the source papers.

Outputs: results/ldsc/h2.tsv, results/ldsc/rg.tsv
"""
import os
import re
import subprocess
from itertools import combinations
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT.parent / "tools"
PY = TOOLS / "ldsc_env" / "bin" / "python"
LDSC = TOOLS / "cbiit_ldsc" / "ldsc.py"
LD = str(ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr") + "/"
MUNGED = ROOT / "data" / "munged"
OUT = ROOT / "results" / "ldsc"

# trait: (sample prevalence, population prevalence); None = report observed scale only
PSYCH = {
    "ASD": (18381 / (18381 + 27969), 0.012),   # Grove 2019, N cases and controls, K from the paper
    "SCZ": (0.5, 0.01),
    "BIP": (0.5, 0.02),
    "MDD": (0.5, 0.15),
    "PTSD": None,                              # z-score meta-analysis mixing symptom scores and diagnoses
    "MDD_EHR": (0.5, 0.15),
    "MDD_Clin": (0.5, 0.15),
}
SUMSTATS = {"ASD": "ASD_Grove2019"}

ENV = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2",
           VECLIB_MAXIMUM_THREADS="2")
NUM = r"(-?[\d.]+(?:e-?\d+)?)"


def path(trait):
    return MUNGED / f"{SUMSTATS.get(trait, trait)}.sumstats.gz"


def ldsc(args, out, done_marker):
    log = Path(str(out) + ".log")
    if log.exists() and done_marker in log.read_text():
        return log.read_text()
    cmd = ["nice", "-n", "10", str(PY), str(LDSC), *args, "--ref-ld-chr", LD, "--w-ld-chr", LD,
           "--out", str(out)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=ENV)
    text = log.read_text()
    if done_marker not in text:
        raise RuntimeError(f"LDSC produced no result for {out.name}; see {log}")
    return text


def grab(text, label):
    m = re.search(re.escape(label) + r":\s*" + NUM + r"(?:\s*\(" + NUM + r"\))?", text)
    return (float(m.group(1)), float(m.group(2)) if m.group(2) else None) if m else (None, None)


def h2_record(trait, text, kind, prev):
    scale = "Liability" if prev else "Observed"
    h2, se = grab(text, f"Total {scale} scale h2")
    icpt, icpt_se = grab(text, "Intercept")
    return {"trait": trait, "kind": kind, "samp_prev": prev[0] if prev else None,
            "pop_prev": prev[1] if prev else None, "scale": scale.lower(), "h2": h2, "h2_se": se,
            "h2_z": h2 / se if h2 is not None and se else None,
            "lambda_gc": grab(text, "Lambda GC")[0], "mean_chi2": grab(text, "Mean Chi^2")[0],
            "intercept": icpt, "intercept_se": icpt_se}


def rg_record(t1, t2, text):
    tail = text[text.find("Summary of Genetic Correlation Results"):].splitlines()
    rec = dict(zip(tail[1].split(), tail[2].split()))
    row = {"trait1": t1, "trait2": t2}
    for k, v in rec.items():
        if k not in ("p1", "p2"):
            row[k] = float(v) if v != "NA" else None
    return row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    psych = [t for t in PSYCH if path(t).exists()]
    endpoints = []
    for done in sorted(MUNGED.glob("*.done")):
        meta = dict(l.split("\t", 1) for l in done.read_text().strip().splitlines())
        if "endpoint" in meta:
            ca, co = int(meta["cases"]), int(meta["controls"])
            endpoints.append((meta["endpoint"], ca / (ca + co)))

    h2_rows = []
    for t in psych:
        prev = PSYCH[t]
        args = ["--h2", str(path(t))] + (["--samp-prev", f"{prev[0]}", "--pop-prev", f"{prev[1]}"] if prev else [])
        text = ldsc(args, OUT / f"h2_{SUMSTATS.get(t, t)}", "scale h2")
        h2_rows.append(h2_record(t, text, "psychiatric", prev))
    for e, p in endpoints:
        text = ldsc(["--h2", str(path(e)), "--samp-prev", f"{p}", "--pop-prev", f"{p}"], OUT / f"h2_{e}",
                    "Total Liability scale h2")
        h2_rows.append(h2_record(e, text, "finngen", (p, p)))

    rg_rows = []
    for t in psych:
        for e, _ in endpoints:
            name = "ASD" if t == "ASD" else t
            text = ldsc(["--rg", f"{path(t)},{path(e)}"], OUT / f"rg_{name}_{e}",
                        "Summary of Genetic Correlation Results")
            rg_rows.append(rg_record(t, e, text))
    for t1, t2 in combinations(psych, 2):
        text = ldsc(["--rg", f"{path(t1)},{path(t2)}"], OUT / f"rg_{t1}_{t2}",
                    "Summary of Genetic Correlation Results")
        rg_rows.append(rg_record(t1, t2, text))

    h2 = pd.DataFrame(h2_rows)
    rg = pd.DataFrame(rg_rows)
    h2.to_csv(OUT / "h2.tsv", sep="\t", index=False)
    rg.to_csv(OUT / "rg.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(h2.round(4).to_string())
    print(rg.round(4).to_string())


if __name__ == "__main__":
    main()
