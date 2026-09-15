"""Stream one FinnGen R12 endpoint, keep HapMap3 SNPs, and write an LDSC-ready file.

The raw FinnGen file (0.3 to 0.8 GB) is never written to disk. It is read straight from
the public bucket, filtered to the HapMap3 rsIDs in the LDSC w_hm3.snplist, and saved as a
small tab-separated file that munge_sumstats.py then turns into .sumstats.gz.

FinnGen R12 columns: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt
af_alt_cases af_alt_controls (GRCh38; beta is for alt). LDSC merges on rsID, so the genome
build does not matter for the regression itself.

Usage: python src/01_fetch_munge_finngen.py ENDPOINT [ENDPOINT ...]
Resumable: an endpoint with a finished marker in data/munged is skipped.
"""
import csv
import gzip
import hashlib
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT.parent / "tools"
PY = TOOLS / "ldsc_env" / "bin" / "python"
MUNGE = TOOLS / "cbiit_ldsc" / "munge_sumstats.py"
HM3 = ROOT / "data" / "ref" / "z8182036" / "eur_w_ld_chr" / "w_hm3.snplist"
MANIFEST = ROOT / "data" / "raw" / "finngen_R12_manifest.tsv"
OUT = ROOT / "data" / "munged"
TMP = ROOT / "data" / "processed" / "hm3_filtered"
LOG = ROOT / "logs"

MAF_MIN = 0.01


def load_manifest():
    with open(MANIFEST) as f:
        return {r["phenocode"]: r for r in csv.DictReader(f, delimiter="\t")}


def load_hm3():
    with open(HM3) as f:
        next(f)
        return {line.split("\t", 1)[0] for line in f}


def stream_filter(url, hm3, dest, keys_out=None):
    """Download and filter in one pass. Returns (rows_read, rows_kept, md5 of the gz stream)."""
    md5 = hashlib.md5()

    class Tee:
        def __init__(self, raw):
            self.raw = raw

        def read(self, n=-1):
            b = self.raw.read(n)
            md5.update(b)
            return b

    n_read = n_kept = 0
    with urllib.request.urlopen(url, timeout=120) as resp, open(dest, "w") as out:
        gz = gzip.GzipFile(fileobj=Tee(resp))
        header = gz.readline().decode().lstrip("#").rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        out.write("SNP\tA1\tA2\tBETA\tSE\tP\tFRQ\n")
        for raw in gz:
            n_read += 1
            f = raw.decode().rstrip("\n").split("\t")
            rsids = f[ix["rsids"]]
            if not rsids:
                continue
            # a position can carry several rsIDs separated by commas
            hit = next((r for r in rsids.split(",") if r in hm3), None)
            if hit is None:
                continue
            af = float(f[ix["af_alt"]])
            if min(af, 1 - af) < MAF_MIN:
                continue
            out.write(f"{hit}\t{f[ix['alt']]}\t{f[ix['ref']]}\t{f[ix['beta']]}\t"
                      f"{f[ix['sebeta']]}\t{f[ix['pval']]}\t{af}\n")
            if keys_out is not None:
                keys_out.write(f"{f[ix['chrom']]}:{f[ix['pos']]}:{f[ix['ref']]}:{f[ix['alt']]}\t{hit}\n")
            n_kept += 1
        # drain so the md5 covers the whole object
        while resp.read(1 << 20):
            pass
    return n_read, n_kept, md5.hexdigest()


def run(endpoint, manifest, hm3):
    done = OUT / f"{endpoint}.done"
    if done.exists():
        print(f"{endpoint}: already done")
        return
    row = manifest[endpoint]
    TMP.mkdir(parents=True, exist_ok=True)
    tmp = TMP / f"{endpoint}.hm3.tsv"
    t0 = time.time()
    # Once, save the FinnGen GRCh38 variant key -> rsID map for HapMap3 SNPs. The FinnGen LD matrix
    # names variants as chr:pos:ref:alt, so this map is what turns it into rsID-keyed LD scores.
    keys_path = ROOT / "data" / "processed" / "finngen_hm3_variant_keys.tsv"
    keys_out = None if keys_path.exists() else open(str(keys_path) + ".partial", "w")
    n_read, n_kept, md5 = stream_filter(row["path_https"], hm3, tmp, keys_out)
    if keys_out is not None:
        keys_out.close()
        Path(str(keys_path) + ".partial").rename(keys_path)
    print(f"{endpoint}: read {n_read:,} variants, kept {n_kept:,} HapMap3 with MAF >= {MAF_MIN}, "
          f"md5 {md5}, {time.time() - t0:.0f}s")
    cmd = [str(PY), str(MUNGE), "--sumstats", str(tmp),
           "--N-cas", row["num_cases"], "--N-con", row["num_controls"],
           "--signed-sumstats", "BETA,0", "--merge-alleles", str(HM3),
           "--out", str(OUT / endpoint)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    tmp.unlink()
    done.write_text(f"endpoint\t{endpoint}\ncases\t{row['num_cases']}\ncontrols\t{row['num_controls']}\n"
                    f"url\t{row['path_https']}\nmd5_gz\t{md5}\nvariants_read\t{n_read}\n"
                    f"hm3_kept\t{n_kept}\nfinished\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def main():
    manifest = load_manifest()
    hm3 = load_hm3()
    OUT.mkdir(parents=True, exist_ok=True)
    for endpoint in sys.argv[1:]:
        run(endpoint, manifest, hm3)


if __name__ == "__main__":
    main()
