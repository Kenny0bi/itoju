## Finnish LD scores for HapMap3 SNPs from the public FinnGen R12 LD matrix.
##
## R implementation of src/06_finnish_ldscores.py. Deterministic: same window, same
## r^2 sum, same self-LD term, so the written .l2.ldscore.gz files match the Python.
##
## The standard LDSC reference (eur_w_ld_chr) sums r^2 over HapMap3 SNPs in 489
## European 1000 Genomes samples. FinnGen is Finnish, a population with its own LD
## structure, so this builds the same quantity from LD measured in the 520,210
## FinnGen R12 samples and lets rg be recomputed with a matched reference.
##
## Source: gs://finngen-public-data-r12/ld_matrix/finngen_r12_chr{c}_ld.tsv.gz, all
## variant pairs within 3 Mb with r^2 > 0.01, listed in both directions
## (PLINK 2 --r-unphased). Each file is streamed and never stored.
##
## For each HapMap3 SNP j:  L2_j = 1 + sum over HapMap3 SNPs k != j within WINDOW_BP of r^2_jk
## The 1 is the SNP's LD with itself, which LDSC includes and the FinnGen file omits.
## Two differences from eur_w_ld_chr, to report: the window is 1 Mb in physical
## distance (no genetic map here) instead of 1 cM, and pairs with r^2 <= 0.01 are
## missing from the source, which biases L2 slightly downward. The comparison against
## eur_w_ld_chr quantifies both together.
##
## Usage: Rscript src/06_finnish_ldscores.R CHR [CHR ...]   (resumable per chromosome)
## Output: data/ref/finngen_ld/{c}.l2.ldscore.gz (CHR SNP BP L2), {c}.l2.M_5_50, and a
## comparison line.

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT <- find_root()
KEYS <- file.path(ROOT, "data", "processed", "finngen_hm3_variant_keys.tsv")
EUR  <- file.path(ROOT, "data", "ref", "z8182036", "eur_w_ld_chr")
OUT  <- file.path(ROOT, "data", "ref", "finngen_ld")
URL  <- "https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/finngen_r12_chr%s_ld.tsv.gz"
WINDOW_BP <- 1000000L

load_keys <- function(chrom) {
  ## HapMap3 keys for one chromosome, rewritten to the LD matrix naming.
  ## Summary statistics use 21:5031188:T:C; the LD matrix uses chr21_5031188_T_C.
  k <- fread(KEYS, sep = "\t", header = FALSE, col.names = c("key", "rsid"), colClasses = "character")
  parts <- tstrsplit(k$key, ":", fixed = TRUE)
  keep <- sub("^chr", "", parts[[1L]]) == as.character(chrom)
  data.table(name = sprintf("chr%s_%s_%s_%s", chrom, parts[[2L]][keep], parts[[3L]][keep], parts[[4L]][keep]),
             rsid = k$rsid[keep], pos = as.integer(parts[[2L]][keep]))
}

build <- function(chrom) {
  done <- file.path(OUT, sprintf("%s.l2.ldscore.gz", chrom))
  if (file.exists(done)) { cat(sprintf("chr%s: already built\n", chrom)); return(invisible(NULL)) }
  keys <- load_keys(chrom)
  setkey(keys, name)
  pos <- setNames(keys$pos, keys$name)
  l2  <- setNames(numeric(length(pos)), names(pos))
  n_lines <- 0L; n_used <- 0L
  t0 <- Sys.time()

  local_gz <- tempfile(fileext = ".gz"); on.exit(unlink(local_gz), add = TRUE)
  utils::download.file(sprintf(URL, chrom), local_gz, mode = "wb", quiet = TRUE)
  con <- gzfile(local_gz, "rt"); on.exit(close(con), add = TRUE)
  header <- strsplit(sub("^#", "", readLines(con, n = 1L)), "\t", fixed = TRUE)[[1L]]
  ix <- setNames(seq_along(header), header)

  repeat {
    chunk <- readLines(con, n = 500000L)
    if (!length(chunk)) break
    n_lines <- n_lines + length(chunk)
    f  <- tstrsplit(chunk, "\t", fixed = TRUE)
    v1 <- f[[ix[["variant1"]]]]; v2 <- f[[ix[["variant2"]]]]
    r2 <- suppressWarnings(as.numeric(f[[ix[["r2"]]]]))
    ok <- v1 %in% names(pos) & v2 %in% names(pos)
    if (!any(ok)) next
    ok[ok] <- abs(pos[v1[ok]] - pos[v2[ok]]) <= WINDOW_BP
    if (!any(ok)) next
    add <- tapply(r2[ok], v1[ok], sum)
    l2[names(add)] <- l2[names(add)] + add
    n_used <- n_used + sum(ok)
  }

  df <- data.table(CHR = as.integer(chrom), SNP = keys$rsid, BP = keys$pos,
                   L2 = 1.0 + as.numeric(l2[keys$name]))
  setorder(df, BP)
  df <- unique(df, by = "SNP")
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

  eur <- fread(file.path(EUR, sprintf("%s.l2.ldscore.gz", chrom)), sep = "\t")
  m <- merge(df, eur[, .(SNP, L2_eur = L2, MAF)], by = "SNP")
  setnames(m, "L2", "L2_fin")
  ## Keep the regression SNP set comparable: only SNPs present in the EUR reference.
  df <- df[SNP %in% eur$SNP]
  fwrite(df, done, sep = "\t", compress = "gzip")

  ## M_5_50: number of SNPs with MAF >= 5% entering the sum. Use the EUR count for
  ## the same SNP set scaled to the SNPs kept here, since per-SNP Finnish MAF is not
  ## carried in the key file.
  m550 <- as.numeric(strsplit(readLines(file.path(EUR, sprintf("%s.l2.M_5_50", chrom)), warn = FALSE)[1L], "\\s+")[[1L]][1L])
  writeLines(sprintf("%.0f", m550 * nrow(df) / nrow(eur)), file.path(OUT, sprintf("%s.l2.M_5_50", chrom)))

  r <- cor(m$L2_fin, m$L2_eur)
  msg <- sprintf("chr%s: %s pair lines, %s HM3-HM3 pairs within 1 Mb, %s SNPs; mean L2 fin %.2f vs eur %.2f; r(fin, eur) = %.4f; %.0fs",
                 chrom, format(n_lines, big.mark = ","), format(n_used, big.mark = ","),
                 format(nrow(df), big.mark = ","), mean(m$L2_fin), mean(m$L2_eur), r,
                 as.numeric(difftime(Sys.time(), t0, units = "secs")))
  writeLines(msg, file.path(OUT, sprintf("%s.compare.txt", chrom)))
  cat(msg, "\n")
}

if (sys.nframe() == 0L) for (c in commandArgs(trailingOnly = TRUE)) build(as.integer(c))
