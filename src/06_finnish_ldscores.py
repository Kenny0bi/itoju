"""Finnish LD scores for HapMap3 SNPs from the public FinnGen R12 LD matrix.

The standard LDSC reference (eur_w_ld_chr) sums r^2 over HapMap3 SNPs in 489 European 1000 Genomes
samples. FinnGen is Finnish, a population with its own LD structure, so this builds the same
quantity from LD measured in the 520,210 FinnGen R12 samples and lets rg be recomputed with a
matched reference.

Source: gs://finngen-public-data-r12/ld_matrix/finngen_r12_chr{c}_ld.tsv.gz, all variant pairs within
3 Mb with r^2 > 0.01, listed in both directions (PLINK 2 --r-unphased). Each file is streamed and
never stored.

For each HapMap3 SNP j:  L2_j = 1 + sum over HapMap3 SNPs k != j within WINDOW_BP of r^2_jk
The 1 is the SNP's LD with itself, which LDSC includes and the FinnGen file omits. Two
differences from eur_w_ld_chr, to report: the window is 1 Mb in physical distance (no genetic map
here) instead of 1 cM, and pairs with r^2 <= 0.01 are missing from the source, which biases L2
slightly downward. The comparison against eur_w_ld_chr quantifies both together.

Usage: python src/06_finnish_ldscores.py CHR [CHR ...]   (resumable per chromosome)
Output: data/ref/finngen_ld/{c}.l2.ldscore.gz (CHR SNP BP L2), {c}.l2.M_5_50, and a comparison line.
"""
import gzip
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEYS = ROOT / "data" / "processed" / "finngen_hm3_variant_keys.tsv"
EUR = ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr"
OUT = ROOT / "data" / "ref" / "finngen_ld"
URL = "https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/finngen_r12_chr{c}_ld.tsv.gz"
WINDOW_BP = 1_000_000


def load_keys(chrom):
    """HapMap3 keys for one chromosome, rewritten to the LD matrix naming.
    Summary statistics use 21:5031188:T:C; the LD matrix uses chr21_5031188_T_C."""
    keys = {}
    with open(KEYS) as f:
        for line in f:
            key, rsid = line.rstrip("\n").split("\t")
            c, pos, ref, alt = key.split(":")
            if c.replace("chr", "") == str(chrom):
                keys[f"chr{chrom}_{pos}_{ref}_{alt}"] = rsid
    return keys


def build(chrom):
    done = OUT / f"{chrom}.l2.ldscore.gz"
    if done.exists():
        print(f"chr{chrom}: already built")
        return
    keys = load_keys(chrom)
    pos = {k: int(k.split("_")[1]) for k in keys}
    l2 = defaultdict(float)
    n_lines = n_used = 0
    t0 = time.time()
    with urllib.request.urlopen(URL.format(c=chrom), timeout=300) as resp:
        gz = gzip.GzipFile(fileobj=resp)
        header = gz.readline().decode().lstrip("#").rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        i1, i2, ir2 = ix["variant1"], ix["variant2"], ix["r2"]
        for raw in gz:
            n_lines += 1
            f = raw.split(b"\t")
            v1 = f[i1].decode()
            if v1 not in keys:
                continue
            v2 = f[i2].decode()
            if v2 not in keys or abs(pos[v1] - pos[v2]) > WINDOW_BP:
                continue
            l2[v1] += float(f[ir2])
            n_used += 1
    rows = [{"CHR": chrom, "SNP": keys[k], "BP": pos[k], "L2": 1.0 + l2.get(k, 0.0)} for k in keys]
    df = pd.DataFrame(rows).sort_values("BP").drop_duplicates("SNP")
    OUT.mkdir(parents=True, exist_ok=True)

    eur = pd.read_csv(EUR / f"{chrom}.l2.ldscore.gz", sep="\t")
    m = df.merge(eur[["SNP", "L2", "MAF"]], on="SNP", suffixes=("_fin", "_eur"))
    # Keep the regression SNP set comparable: only SNPs present in the EUR reference.
    df = df[df.SNP.isin(eur.SNP)]
    df.to_csv(done, sep="\t", index=False, compression="gzip")
    # M_5_50: number of SNPs with MAF >= 5% entering the sum. Use the EUR count for the same SNP set
    # scaled to the SNPs kept here, since per-SNP Finnish MAF is not carried in the key file.
    m550 = float(open(EUR / f"{chrom}.l2.M_5_50").read().split()[0])
    (OUT / f"{chrom}.l2.M_5_50").write_text(f"{m550 * len(df) / len(eur):.0f}\n")
    r = np.corrcoef(m.L2_fin, m.L2_eur)[0, 1]
    msg = (f"chr{chrom}: {n_lines:,} pair lines, {n_used:,} HM3-HM3 pairs within 1 Mb, {len(df):,} SNPs; "
           f"mean L2 fin {m.L2_fin.mean():.2f} vs eur {m.L2_eur.mean():.2f}; r(fin, eur) = {r:.4f}; "
           f"{time.time() - t0:.0f}s")
    (OUT / f"{chrom}.compare.txt").write_text(msg + "\n")
    print(msg, flush=True)


if __name__ == "__main__":
    for c in sys.argv[1:]:
        build(int(c))
