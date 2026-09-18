## Prepare the PGC psychiatric GWAS for LDSC and check allele alignment against 1000 Genomes.
##
## R implementation of src/04_prep_pgc.py. Deterministic throughout, so the written
## tables and the reported allele-frequency correlations match the Python exactly.
## munge_sumstats.py has no R equivalent and is invoked identically by both.
##
## Each file is streamed (never modified; the originals are read-only symlinks), the
## '##' metadata block is skipped, rows are restricted to HapMap3 rsIDs, and a clean
## table is written with SNP A1 A2 SIGNED P N FRQ [INFO]. munge_sumstats.py then
## produces data/munged/<trait>.sumstats.gz.
##
## Alignment check. LDSC reads the signed statistic as the effect of A1. A mislabelled
## allele column flips the sign of every z-score, which flips the sign of every rg. The
## reported allele frequency (FCON or FREQ) belongs to the same allele as the effect, so
## it should track the 1000 Genomes EUR frequency of that allele closely (r near +1). A
## swap would show r near -1.
##
## Usage: Rscript src/04_prep_pgc.R TRAIT [TRAIT ...]   (resumable via data/munged/<trait>.done)

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT  <- find_root()
TOOLS <- normalizePath(file.path(ROOT, "..", "tools"), mustWork = FALSE)
PY    <- file.path(TOOLS, "ldsc_env", "bin", "python")
MUNGE <- file.path(TOOLS, "cbiit_ldsc", "munge_sumstats.py")
HM3   <- file.path(ROOT, "data", "ref", "z8182036", "eur_w_ld_chr", "w_hm3.snplist")
FRQ   <- file.path(ROOT, "data", "ref", "1000G_Phase3_frq")
RAW   <- file.path(ROOT, "data", "raw")
TMP   <- file.path(ROOT, "data", "processed", "pgc")
OUT   <- file.path(ROOT, "data", "munged")

## n: column name, or a function of the row
TRAITS <- list(
  SCZ = list(file = "pgc/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
             snp = "ID", a1 = "A1", a2 = "A2", signed = "BETA", p = "PVAL", n = "NEFF",
             frq = "FCON", info = "IMPINFO"),
  BIP = list(file = "pgc/pgc-bip2021-all.vcf.tsv.gz",
             snp = "ID", a1 = "A1", a2 = "A2", signed = "BETA", p = "PVAL",
             n = function(d) 2 * as.numeric(d[["NEFFDIV2"]]), frq = "FCON", info = "IMPINFO"),
  MDD = list(file = "pgc/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz",
             snp = "ID", a1 = "EA", a2 = "NEA", signed = "BETA", p = "PVAL", n = "NEFF",
             frq = "FCON", info = "IMPINFO"),
  PTSD = list(file = "pgc/eur_ptsd_pcs_v4_aug3_2021.vcf.gz",
              snp = "ID", a1 = "A1", a2 = "A2", signed = "Z", p = "P", n = "NEFF",
              frq = "FREQ", info = NULL),
  MDD_EHR = list(file = "mdd2025/pgc-mdd2025_EHR_eur_v3-49-24-11.tsv.gz",
                 snp = "ID", a1 = "EA", a2 = "NEA", signed = "BETA", p = "PVAL", n = "NEFF",
                 frq = "FCON", info = "IMPINFO"),
  MDD_Clin = list(file = "mdd2025/pgc-mdd2025_Clin_eur_v3-49-24-11.tsv.gz",
                  snp = "ID", a1 = "EA", a2 = "NEA", signed = "BETA", p = "PVAL", n = "NEFF",
                  frq = "FCON", info = "IMPINFO")
)
COMP <- c(A = "T", T = "A", C = "G", G = "C")

load_hm3 <- function() unique(fread(HM3, sep = "\t", select = 1L, colClasses = "character")[[1L]])

load_frq <- function(hm3) {
  ## rsID -> (A1, A2, frequency of A1) in 1000 Genomes phase 3 EUR, HapMap3 SNPs only.
  ## Restricting to HapMap3 keeps the table near 1.2M entries instead of the
  ## genome-wide set.
  parts <- lapply(1:22, function(c) {
    d <- fread(file.path(FRQ, sprintf("1000G.EUR.QC.%d.frq", c)))
    setnames(d, c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS")[seq_len(ncol(d))])
    d[SNP %in% hm3, .(SNP, A1, A2, MAF)]
  })
  ref <- rbindlist(parts)
  setkey(ref, SNP)
  ref[]
}

count_meta_lines <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  n <- 0L
  repeat {
    line <- readLines(con, n = 1L)
    if (!length(line) || !startsWith(line, "##")) break
    n <- n + 1L
  }
  n
}

prep <- function(trait, spec, hm3, ref) {
  src <- file.path(RAW, spec$file)
  dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(TMP, paste0(trait, ".hm3.tsv"))

  d <- fread(src, sep = "\t", skip = count_meta_lines(src))
  n_read <- nrow(d)
  d <- d[get(spec$snp) %in% hm3]
  a1 <- toupper(d[[spec$a1]]); a2 <- toupper(d[[spec$a2]])
  nvals <- if (is.function(spec$n)) spec$n(d) else d[[spec$n]]

  outdt <- data.table(SNP = d[[spec$snp]], A1 = a1, A2 = a2,
                      SIGNED = d[[spec$signed]], P = d[[spec$p]], N = nvals, FRQ = d[[spec$frq]])
  if (!is.null(spec$info)) outdt[, INFO := d[[spec$info]]]
  fwrite(outdt, tmp, sep = "\t")
  n_kept <- nrow(outdt)

  ## palindromic SNPs cannot be checked by allele label
  keep <- COMP[a1] != a2
  keep[is.na(keep)] <- FALSE
  m <- match(outdt$SNP, ref$SNP)
  ok <- keep & !is.na(m)
  gw_f <- numeric(0); ref_f <- numeric(0)
  if (any(ok)) {
    idx <- which(ok)
    r1 <- ref$A1[m[idx]]; r2 <- ref$A2[m[idx]]; rmaf <- ref$MAF[m[idx]]
    fa <- suppressWarnings(as.numeric(outdt$FRQ[idx]))
    direct <- a1[idx] == r1 & a2[idx] == r2
    swapd  <- a1[idx] == r2 & a2[idx] == r1
    gw_f  <- c(fa[direct], fa[swapd])
    ref_f <- c(rmaf[direct], 1 - rmaf[swapd])
    good  <- is.finite(gw_f) & is.finite(ref_f)
    gw_f <- gw_f[good]; ref_f <- ref_f[good]
  }
  r_frq <- if (length(gw_f) > 2L) cor(gw_f, ref_f) else NA_real_
  far   <- if (length(gw_f)) mean(abs(gw_f - ref_f) > 0.2) else NA_real_
  list(tmp = tmp, n_read = n_read, n_kept = n_kept, n_chk = length(gw_f), r_frq = r_frq, far = far)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  hm3 <- load_hm3()
  ref <- load_frq(hm3)
  for (trait in args) {
    spec <- TRAITS[[trait]]
    done <- file.path(OUT, paste0(trait, ".done"))
    if (file.exists(done)) { cat(trait, ": already done\n", sep = ""); next }
    if (!file.exists(file.path(RAW, spec$file))) { cat(trait, ": source missing, skipped\n", sep = ""); next }
    t0 <- Sys.time()
    r  <- prep(trait, spec, hm3, ref)
    cat(sprintf("%s: read %s, kept %s HM3; allele-frequency check on %s non-palindromic SNPs: r = %.4f, |diff| > 0.2 in %.4f%%; %.0fs\n",
                trait, format(r$n_read, big.mark = ","), format(r$n_kept, big.mark = ","),
                format(r$n_chk, big.mark = ","), r$r_frq, 100 * r$far,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (!is.na(r$r_frq) && r$r_frq < 0.9)
      stop(sprintf("%s: allele frequencies do not track 1000G EUR (r = %.3f); stop and inspect", trait, r$r_frq),
           call. = FALSE)
    cmd <- c(MUNGE, "--sumstats", r$tmp, "--signed-sumstats", "SIGNED,0",
             "--N-col", "N", "--frq", "FRQ", "--merge-alleles", HM3, "--out", file.path(OUT, trait))
    cmd <- if (!is.null(spec$info)) c(cmd, "--info", "INFO", "--info-min", "0.9") else c(cmd, "--ignore", "INFO")
    system2("nice", c("-n", "10", PY, cmd), stdout = FALSE)
    unlink(r$tmp)
    writeLines(c(paste0("trait\t", trait), paste0("source\t", spec$file),
                 paste0("hm3_kept\t", r$n_kept), paste0("frq_check_n\t", r$n_chk),
                 sprintf("frq_check_r\t%.6f", r$r_frq), sprintf("frq_far_frac\t%.6f", r$far),
                 paste0("finished\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), done)
  }
}

if (sys.nframe() == 0L) main()
