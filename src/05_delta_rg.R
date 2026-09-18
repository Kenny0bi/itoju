## Standard error of a difference between two genetic correlations that share a trait.
## The two estimates share people, so their errors are correlated and subtracting
## block by block cancels the shared noise.
##
## Usage: Rscript src/05_delta_rg.R LOG_A LOG_B

Z_MDD <- qnorm(0.975) + qnorm(0.80)

load_pair <- function(log_path) {
  text <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  tail_lines <- strsplit(substring(text, regexpr("Summary of Genetic Correlation Results", text)), "\n")[[1L]]
  hdr <- strsplit(trimws(tail_lines[2L]), "\\s+")[[1L]]
  val <- strsplit(trimws(tail_lines[3L]), "\\s+")[[1L]]
  rec <- setNames(as.list(val), hdr)

  stem <- sub("\\.log$", "", basename(log_path))
  dir  <- dirname(log_path)
  deletes <- list()
  for (part in c("hsq1", "hsq2", "gencov")) {
    hits <- list.files(dir, pattern = paste0("^", stem, ".*\\.", part, "\\.delete$"), full.names = TRUE)
    if (length(hits) != 1L)
      stop(sprintf("expected one %s.delete for %s, found %d", part, basename(log_path), length(hits)),
           call. = FALSE)
    deletes[[part]] <- as.numeric(readLines(hits[1L], warn = FALSE))
  }
  c(list(rg = as.numeric(rec$rg), se = as.numeric(rec$se), p2 = rec$p2), deletes)
}

pseudovalues <- function(pair) {
  n      <- length(pair$gencov)
  rg_del <- pair$gencov / sqrt(pair$hsq1 * pair$hsq2)
  n * pair$rg - (n - 1) * rg_del
}

jackknife_se <- function(pseudo) sqrt(var(pseudo) / length(pseudo))

compare <- function(log_a, log_b) {
  a <- load_pair(log_a); b <- load_pair(log_b)
  if (length(a$gencov) != length(b$gencov))
    stop("block counts differ; the pairs were not run on the same SNP set", call. = FALSE)
  pa <- pseudovalues(a); pb <- pseudovalues(b)
  for (nm in c("A", "B")) {
    pair   <- if (nm == "A") a else b
    pseudo <- if (nm == "A") pa else pb
    rebuilt <- jackknife_se(pseudo)
    gap     <- abs(rebuilt - pair$se) / pair$se
    ## LDSC prints SE to 4 decimals, so allow for rounding in the printed value
    if (abs(rebuilt - pair$se) > 6e-5 && gap > 1e-3)
      stop(sprintf("pair %s: rebuilt SE %.5f does not match LDSC SE %.4f", nm, rebuilt, pair$se),
           call. = FALSE)
  }
  delta    <- a$rg - b$rg
  se_delta <- jackknife_se(pa - pb)
  list(rg_a = a$rg, se_a = a$se, rg_b = b$rg, se_b = b$se,
       delta = delta, se_delta = se_delta, z = delta / se_delta,
       p = 2 * pnorm(abs(delta / se_delta), lower.tail = FALSE),
       error_corr = cor(pa, pb),
       se_delta_if_independent = sqrt(a$se^2 + b$se^2),
       min_detectable_delta = Z_MDD * se_delta,
       n_blocks = length(pa),
       rebuilt_se_a = jackknife_se(pa), rebuilt_se_b = jackknife_se(pb))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) stop("usage: Rscript src/05_delta_rg.R LOG_A LOG_B", call. = FALSE)
  out <- compare(args[1L], args[2L])
  for (k in names(out)) {
    v <- out[[k]]
    cat(sprintf(if (is.numeric(v) && !is.integer(v)) "%-26s %.5f\n" else "%-26s %s\n", k, v))
  }
}
