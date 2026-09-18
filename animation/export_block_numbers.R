## Real numbers for the jackknife part of the animation.
##
## R implementation of animation/export_block_numbers.py. Deterministic: it reads
## the same LDSC .delete files and the same pairwise_delta.tsv, so it writes the
## same block_numbers.json to floating-point precision.
##
## LDSC's standard error comes from a jackknife over 200 genome blocks. On the common
## SNP set, block k is the same stretch of genome for every trait, so for two
## definitions of one condition the per-block pseudovalues can be paired. This
## exports, for depression with sleep apnoea (hospital records) and with any sleep
## disorder, the 200 pseudovalues under each definition, both estimates, the paired
## standard error of the difference, and the standard error if the pairing is broken
## (as if the two estimates shared no people). Values come from
## results/definition_effect/ldsc, the same files the paper uses.
##
## Writes animation/block_numbers.json.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT <- find_root()
RUNS <- file.path(ROOT, "results", "definition_effect", "ldsc")

pseudo <- function(first, second, tag = "") {
  log_path <- file.path(RUNS, sprintf("rg_%s%s.log", first, tag))
  text  <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  lines <- strsplit(substring(text, regexpr("Summary of Genetic Correlation Results", text)), "\n")[[1L]][-1L]
  hdr   <- strsplit(trimws(lines[1L]), "\\s+")[[1L]]
  rec   <- NULL
  for (l in lines[-1L]) {
    v <- strsplit(trimws(l), "\\s+")[[1L]]
    if (length(v) != length(hdr)) next
    r <- setNames(as.list(v), hdr)
    if (!is.null(r$p2) && endsWith(r$p2, sprintf("/%s.sumstats.gz", second))) { rec <- r; break }
  }
  if (is.null(rec)) stop(sprintf("no row for %s in %s", second, basename(log_path)), call. = FALSE)

  stem <- sprintf("%s%s.sumstats.gz_%s.sumstats.gz", tools::file_path_sans_ext(basename(log_path)), first, second)
  g  <- as.numeric(readLines(file.path(RUNS, paste0(stem, ".gencov.delete")), warn = FALSE))
  h1 <- as.numeric(readLines(file.path(RUNS, paste0(stem, ".hsq1.delete")),   warn = FALSE))
  h2 <- as.numeric(readLines(file.path(RUNS, paste0(stem, ".hsq2.delete")),   warn = FALSE))
  n  <- length(g)
  rg <- as.numeric(rec$rg)
  list(pv = n * rg - (n - 1) * g / sqrt(h1 * h2), rg = rg, se = as.numeric(rec$se))
}

main <- function() {
  a <- pseudo("MDD", "G6_SLEEPAPNO")
  b <- pseudo("MDD", "SLEEP")
  n <- length(a$pv)
  se_pair   <- sqrt(var(b$pv - a$pv) / n)
  se_broken <- sqrt((var(a$pv) + var(b$pv)) / n)

  pw  <- fread(file.path(ROOT, "results", "definition_effect", "pairwise_delta.tsv"), sep = "\t")
  row <- pw[shared == "MDD" & def1 == "G6_SLEEPAPNO" & def2 == "SLEEP"][1L]
  stopifnot(abs(se_pair - row$se_delta) < 1e-6)
  stopifnot(abs((b$rg - a$rg) + row$delta) < 1e-3)

  out <- list(trait = "MDD", def1 = "G6_SLEEPAPNO", def2 = "SLEEP", n_blocks = n,
              rg1 = a$rg, rg2 = b$rg, se1 = a$se, se2 = b$se,
              shift = b$rg - a$rg, se_paired = se_pair, se_broken = se_broken,
              p = row$p, error_corr = row$error_corr, min_detectable = row$min_detectable_delta,
              dev1 = round(a$pv - a$rg, 4), dev2 = round(b$pv - b$rg, 4))

  ## which chromosome each jackknife block starts on: blocks are contiguous runs of
  ## the common SNP set in genome order
  snps <- fread(file.path(ROOT, "data", "munged_common", "MDD.sumstats.gz"), sep = "\t", select = "SNP")
  ld <- rbindlist(lapply(1:22, function(c)
    fread(file.path(ROOT, "data", "ref", "z8182036", "eur_w_ld_chr", sprintf("%d.l2.ldscore.gz", c)),
          sep = "\t", select = c("CHR", "SNP", "BP"))))
  order_dt <- merge(snps, ld, by = "SNP")
  setorder(order_dt, CHR, BP)
  m <- nrow(order_dt)
  starts <- ((seq_len(n) - 1L) * m) %/% n + 1L
  out$n_snps_common <- m
  out$block_chrom   <- as.integer(order_dt$CHR[starts])

  write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = TRUE),
        file.path(ROOT, "animation", "block_numbers.json"))
  cat(sprintf("blocks %d; rg %.4f -> %.4f, shift %+.4f; SE paired %.4f, broken %.4f; error corr %.3f; p %.2e; pseudovalue sd %.3f, %.3f\n",
              n, a$rg, b$rg, b$rg - a$rg, se_pair, se_broken, row$error_corr, row$p, sd(a$pv), sd(b$pv)))
}

if (sys.nframe() == 0L) main()
