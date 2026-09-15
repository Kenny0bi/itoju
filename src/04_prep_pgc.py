"""Prepare the PGC psychiatric GWAS for LDSC and check allele alignment against 1000 Genomes.

Each file is streamed (never modified; the originals are read-only symlinks), the '##' metadata
block is skipped, rows are restricted to HapMap3 rsIDs, and a clean table is written with
SNP A1 A2 SIGNED P N FRQ [INFO]. munge_sumstats.py then produces data/munged/<trait>.sumstats.gz.

Alignment check. LDSC reads the signed statistic as the effect of A1. A mislabelled allele column
flips the sign of every z-score, which flips the sign of every rg. The reported allele frequency
(FCON or FREQ) belongs to the same allele as the effect, so it should track the 1000 Genomes EUR
frequency of that allele closely (r near +1). A swap would show r near -1.

Usage: python src/04_prep_pgc.py TRAIT [TRAIT ...]   (resumable via data/munged/<trait>.done)
"""
import gzip
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT.parent / "tools"
PY = TOOLS / "ldsc_env" / "bin" / "python"
MUNGE = TOOLS / "cbiit_ldsc" / "munge_sumstats.py"
HM3 = ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr" / "w_hm3.snplist"
FRQ = ROOT / "data" / "ref" / "1000G_Phase3_frq"
RAW = ROOT / "data" / "raw"
TMP = ROOT / "data" / "processed" / "pgc"
OUT = ROOT / "data" / "munged"

# n: column name, or a function of the row dict
TRAITS = {
    "SCZ": dict(file="pgc/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", snp="ID", a1="A1", a2="A2",
                signed="BETA", p="PVAL", n="NEFF", frq="FCON", info="IMPINFO"),
    "BIP": dict(file="pgc/pgc-bip2021-all.vcf.tsv.gz", snp="ID", a1="A1", a2="A2",
                signed="BETA", p="PVAL", n=lambda r: 2 * float(r["NEFFDIV2"]), frq="FCON", info="IMPINFO"),
    "MDD": dict(file="pgc/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz", snp="ID", a1="EA", a2="NEA",
                signed="BETA", p="PVAL", n="NEFF", frq="FCON", info="IMPINFO"),
    "PTSD": dict(file="pgc/eur_ptsd_pcs_v4_aug3_2021.vcf.gz", snp="ID", a1="A1", a2="A2",
                 signed="Z", p="P", n="NEFF", frq="FREQ", info=None),
    "MDD_EHR": dict(file="mdd2025/pgc-mdd2025_EHR_eur_v3-49-24-11.tsv.gz", snp="ID", a1="EA", a2="NEA",
                    signed="BETA", p="PVAL", n="NEFF", frq="FCON", info="IMPINFO"),
    "MDD_Clin": dict(file="mdd2025/pgc-mdd2025_Clin_eur_v3-49-24-11.tsv.gz", snp="ID", a1="EA", a2="NEA",
                     signed="BETA", p="PVAL", n="NEFF", frq="FCON", info="IMPINFO"),
}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_hm3():
    with open(HM3) as f:
        next(f)
        return {l.split("\t", 1)[0] for l in f}


def load_frq(hm3):
    """rsID -> (A1, A2, frequency of A1) in 1000 Genomes phase 3 EUR, HapMap3 SNPs only.
    Restricting to HapMap3 keeps the table near 1.2M entries instead of the genome-wide set."""
    ref = {}
    for c in range(1, 23):
        with open(FRQ / f"1000G.EUR.QC.{c}.frq") as f:
            next(f)
            for line in f:
                _, snp, a1, a2, maf, _ = line.split()
                if snp in hm3:
                    ref[snp] = (a1, a2, float(maf))
    return ref


def prep(trait, spec, hm3, ref):
    src = RAW / spec["file"]
    TMP.mkdir(parents=True, exist_ok=True)
    tmp = TMP / f"{trait}.hm3.tsv"
    n_read = n_kept = 0
    gw_f, ref_f = [], []
    with gzip.open(src, "rt") as fh, open(tmp, "w") as out:
        for line in fh:
            if not line.startswith("##"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        ix = {c: i for i, c in enumerate(header)}
        cols = ["SNP", "A1", "A2", "SIGNED", "P", "N", "FRQ"] + (["INFO"] if spec["info"] else [])
        out.write("\t".join(cols) + "\n")
        for line in fh:
            n_read += 1
            f = line.rstrip("\n").split("\t")
            snp = f[ix[spec["snp"]]]
            if snp not in hm3:
                continue
            row = {c: f[i] for c, i in ix.items()}
            a1, a2 = row[spec["a1"]].upper(), row[spec["a2"]].upper()
            n = spec["n"](row) if callable(spec["n"]) else row[spec["n"]]
            vals = [snp, a1, a2, row[spec["signed"]], row[spec["p"]], str(n), row[spec["frq"]]]
            if spec["info"]:
                vals.append(row[spec["info"]])
            out.write("\t".join(vals) + "\n")
            n_kept += 1
            r = ref.get(snp)
            if r and COMP.get(a1) != a2:  # palindromic SNPs cannot be checked by allele label
                if (a1, a2) == (r[0], r[1]):
                    gw_f.append(float(row[spec["frq"]])); ref_f.append(r[2])
                elif (a1, a2) == (r[1], r[0]):
                    gw_f.append(float(row[spec["frq"]])); ref_f.append(1 - r[2])
    gw_f, ref_f = np.array(gw_f), np.array(ref_f)
    r_frq = float(np.corrcoef(gw_f, ref_f)[0, 1])
    far = float(np.mean(np.abs(gw_f - ref_f) > 0.2))
    return tmp, n_read, n_kept, len(gw_f), r_frq, far


def main():
    hm3 = load_hm3()
    ref = load_frq(hm3)
    for trait in sys.argv[1:]:
        spec = TRAITS[trait]
        done = OUT / f"{trait}.done"
        if done.exists():
            print(f"{trait}: already done")
            continue
        if not (RAW / spec["file"]).exists():
            print(f"{trait}: source missing, skipped")
            continue
        t0 = time.time()
        tmp, n_read, n_kept, n_chk, r_frq, far = prep(trait, spec, hm3, ref)
        print(f"{trait}: read {n_read:,}, kept {n_kept:,} HM3; allele-frequency check on {n_chk:,} "
              f"non-palindromic SNPs: r = {r_frq:.4f}, |diff| > 0.2 in {far:.4%}; {time.time() - t0:.0f}s",
              flush=True)
        if r_frq < 0.9:
            raise SystemExit(f"{trait}: allele frequencies do not track 1000G EUR (r = {r_frq:.3f}); stop and inspect")
        cmd = [str(PY), str(MUNGE), "--sumstats", str(tmp), "--signed-sumstats", "SIGNED,0",
               "--N-col", "N", "--frq", "FRQ", "--merge-alleles", str(HM3), "--out", str(OUT / trait)]
        if spec["info"]:
            cmd += ["--info", "INFO", "--info-min", "0.9"]
        else:
            cmd += ["--ignore", "INFO"]
        subprocess.run(["nice", "-n", "10", *cmd], check=True, stdout=subprocess.DEVNULL)
        tmp.unlink()
        done.write_text(f"trait\t{trait}\nsource\t{spec['file']}\nhm3_kept\t{n_kept}\nfrq_check_n\t{n_chk}\n"
                        f"frq_check_r\t{r_frq:.6f}\nfrq_far_frac\t{far:.6f}\n"
                        f"finished\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")


if __name__ == "__main__":
    main()
