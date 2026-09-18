## The core test: does the genetic correlation move when only the case definition changes?
##
## R implementation of src/07_definition_effect.py. LDSC is driven by the same
## subprocess calls and the same .delete files are read back, so the pseudovalues,
## deltas, heterogeneity statistics and SE rebuild checks all match the Python.
##
## Design, fixed in docs/ANALYSIS_PLAN.md before any rg was estimated:
##   DEFINITION_SETS  alternative definitions of one construct (they overlap or nest)
##   NEIGHBOURS       different conditions, kept as a comparison for how far rg moves
##                    between conditions
##   SHARED           the trait each definition is correlated with (autism is focal)
##   POSITIVE CONTROL MDD_EHR vs MDD_Clin, two definitions of depression with a known
##                    definition effect
##
## Steps
##   1. Restrict every trait to one common SNP set (non-missing Z in all), so LDSC's
##      200 jackknife blocks are identical across pairs.
##   2. For each shared trait X, run one LDSC call X,D1,D2,... with --two-step 99999
##      and --print-delete-vals (one-step estimator, evenly spaced blocks).
##   3. Rebuild per-block rg pseudovalues exactly as LDSC's RatioJackknife does, and
##      check that they reproduce LDSC's printed SE.
##   4. Within each definition set and shared trait:
##        - pairwise delta rg, SE from the pseudovalue difference, z, p, error
##          correlation, minimum detectable delta (alpha 0.05, power 0.80);
##        - heterogeneity Q = r' W r - (1' W r)^2 / (1' W 1), W = inverse jackknife
##          covariance, df = k - 1 (pseudo-inverse if the covariance is near singular;
##          condition number reported).
##   5. Join pairwise deltas to the definition distances from src/02_definition_anatomy.
##
## Outputs (results/definition_effect/): rg_common.tsv, pairwise_delta.tsv,
## heterogeneity.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
})

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT   <- find_root()
TOOLS  <- normalizePath(file.path(ROOT, "..", "tools"), mustWork = FALSE)
PY     <- file.path(TOOLS, "ldsc_env", "bin", "python")
LDSC   <- file.path(TOOLS, "cbiit_ldsc", "ldsc.py")
LD     <- paste0(file.path(ROOT, "data", "ref", "z8182036", "eur_w_ld_chr"), "/")
MUNGED <- file.path(ROOT, "data", "munged")
COMMON <- file.path(ROOT, "data", "munged_common")
OUT    <- file.path(ROOT, "results", "definition_effect")
RUNS   <- file.path(OUT, "ldsc")

DEFINITION_SETS <- list(
  Epilepsy        = c("G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"),
  `Sleep apnoea`  = c("G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"),
  Insomnia        = c("F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"),
  Constipation    = c("K11_CONSTIPATION", "K11_OTHFUNC"),
  ADHD            = c("F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"),
  `Intellectual disability` = c("F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE")
)
NEIGHBOURS <- c("G6_STATUSEPI", "G6_SLEEPDISOTH", "F5_SLEEP_NOS", "K11_IBS", "K11_FUNCDYSP",
                "K11_REFLUX", "N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD",
                "KRA_PSY_DEVWIDE_EXMORE", "KRA_PSY_AUTISM_EXMORE")
SHARED  <- c("ASD", "SCZ", "BIP", "MDD", "PTSD")
CONTROL <- c("MDD_EHR", "MDD_Clin")
FILES   <- c(ASD = "ASD_Grove2019")

ENV   <- c("OMP_NUM_THREADS=2", "OPENBLAS_NUM_THREADS=2", "MKL_NUM_THREADS=2", "VECLIB_MAXIMUM_THREADS=2")
Z_MDD <- qnorm(0.975) + qnorm(0.80)

src_path <- function(t) file.path(MUNGED, paste0(if (t %in% names(FILES)) FILES[[t]] else t, ".sumstats.gz"))

build_common <- function(traits) {
  ## Write copies of every trait restricted to SNPs with non-missing Z in all of them.
  dir.create(COMMON, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(COMMON, "snps.txt")
  frames <- lapply(traits, function(t) { d <- fread(src_path(t), sep = "\t"); d[!is.na(Z)] })
  names(frames) <- traits
  common <- Reduce(intersect, lapply(frames, function(f) f$SNP))
  if (file.exists(marker) && setequal(readLines(marker, warn = FALSE), common)) return(length(common))
  for (t in traits) fwrite(frames[[t]][SNP %in% common], file.path(COMMON, paste0(t, ".sumstats.gz")),
                           sep = "\t", compress = "gzip")
  writeLines(sort(common), marker)
  length(common)
}

run_rg <- function(first, others, tag = "") {
  ## One LDSC call: first against each of others, one-step, delete values printed.
  ## `tag` keeps runs with the same first trait but different partners from sharing
  ## an output name (the positive-control runs would otherwise collide).
  dir.create(RUNS, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(RUNS, sprintf("rg_%s%s", first, tag))
  log_path <- paste0(out, ".log")
  if (file.exists(log_path)) {
    txt <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
    n_done <- lengths(regmatches(txt, gregexpr("Genetic Correlation: ", txt, fixed = TRUE)))
    if (n_done >= length(others) && grepl("Summary of Genetic Correlation Results", txt, fixed = TRUE))
      return(log_path)
  }
  paths <- paste(vapply(c(first, others), function(t) file.path(COMMON, paste0(t, ".sumstats.gz")),
                        character(1L)), collapse = ",")
  system2("nice", c("-n", "10", PY, LDSC, "--rg", paths, "--ref-ld-chr", LD, "--w-ld-chr", LD,
                    "--two-step", "99999", "--print-delete-vals", "--out", out),
          stdout = FALSE, stderr = FALSE, env = ENV)
  log_path
}

parse_run <- function(first, others, log_path) {
  ## Summary table plus pseudovalues for each pair in one LDSC run.
  text  <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  lines <- strsplit(substring(text, regexpr("Summary of Genetic Correlation Results", text)), "\n")[[1L]][-1L]
  hdr   <- strsplit(trimws(lines[1L]), "\\s+")[[1L]]
  rows  <- list(); pseudo <- list()
  stem_base <- tools::file_path_sans_ext(basename(log_path))
  for (i in seq_along(others)) {
    t   <- others[i]
    v   <- strsplit(trimws(lines[1L + i]), "\\s+")[[1L]]
    rec <- setNames(as.list(v), hdr)
    rg  <- if (rec$rg == "NA") NA_real_ else as.numeric(rec$rg)
    se  <- if (rec$se == "NA") NA_real_ else as.numeric(rec$se)
    row <- list(trait1 = first, trait2 = t, rg = rg, se = se,
                p = if (rec$p == "NA") NA_real_ else as.numeric(rec$p),
                h2_obs_2 = as.numeric(rec$h2_obs), gcov_int = as.numeric(rec$gcov_int),
                gcov_int_se = as.numeric(rec$gcov_int_se),
                se_rebuilt = NA_real_, se_check_ok = NA)
    if (is.finite(rg)) {
      stem <- sprintf("%s%s.sumstats.gz_%s.sumstats.gz", stem_base, first, t)
      g  <- as.numeric(readLines(file.path(RUNS, paste0(stem, ".gencov.delete")), warn = FALSE))
      h1 <- as.numeric(readLines(file.path(RUNS, paste0(stem, ".hsq1.delete")),   warn = FALSE))
      h2 <- as.numeric(readLines(file.path(RUNS, paste0(stem, ".hsq2.delete")),   warn = FALSE))
      n  <- length(g)
      pv <- n * rg - (n - 1) * g / sqrt(h1 * h2)
      rebuilt <- sqrt(var(pv) / n)
      row$se_rebuilt  <- rebuilt
      row$se_check_ok <- is.finite(rebuilt) && (abs(rebuilt - se) <= 6e-5 || abs(rebuilt - se) / se < 1e-3)
      pseudo[[t]] <- pv
    }
    rows[[length(rows) + 1L]] <- as.data.table(row)
  }
  list(rows = rbindlist(rows), pseudo = pseudo)
}

pairwise <- function(x, family, defs, rg, pseudo) {
  out <- list()
  defs <- defs[defs %in% names(pseudo)]
  if (length(defs) < 2L) return(out)
  cb <- combn(defs, 2L)
  for (i in seq_len(ncol(cb))) {
    a <- cb[1L, i]; b <- cb[2L, i]
    d    <- rg[[a]] - rg[[b]]
    diff <- pseudo[[a]] - pseudo[[b]]
    n    <- length(diff)
    se   <- sqrt(var(diff) / n)
    out[[length(out) + 1L]] <- data.table(
      shared = x, family = family, def1 = a, def2 = b, rg1 = rg[[a]], rg2 = rg[[b]],
      delta = d, se_delta = se, z = d / se, p = 2 * pnorm(abs(d / se), lower.tail = FALSE),
      error_corr = cor(pseudo[[a]], pseudo[[b]]), min_detectable_delta = Z_MDD * se)
  }
  out
}

heterogeneity <- function(x, family, defs, rg, pseudo) {
  defs <- defs[defs %in% names(pseudo)]
  if (length(defs) < 2L) return(NULL)
  r <- vapply(defs, function(d) rg[[d]], numeric(1L))
  P <- do.call(cbind, pseudo[defs])
  S <- cov(P) / nrow(P)
  cond <- kappa(S, exact = TRUE)
  W <- MASS::ginv(S)
  one <- rep(1, length(r))
  pooled <- as.numeric(one %*% W %*% r / (one %*% W %*% one))
  Q <- as.numeric(r %*% W %*% r - (one %*% W %*% r)^2 / (one %*% W %*% one))
  df <- length(r) - 1L
  data.table(shared = x, family = family, k = length(r), pooled_rg = pooled, Q = Q, df = df,
             p = pchisq(Q, df, lower.tail = FALSE), rg_range = max(r) - min(r),
             cov_condition_number = cond)
}

main <- function() {
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  endpoints <- c(unlist(DEFINITION_SETS, use.names = FALSE), NEIGHBOURS)
  traits <- c(Filter(function(t) file.exists(src_path(t)), c(SHARED, CONTROL)),
              Filter(function(e) file.exists(src_path(e)), endpoints))
  missing <- Filter(function(t) !file.exists(src_path(t)), c(SHARED, CONTROL, endpoints))
  if (length(missing)) cat("not yet munged, skipped:", paste(missing, collapse = " "), "\n")
  n_common <- build_common(traits)
  cat(sprintf("common SNP set: %s SNPs across %d traits\n", format(n_common, big.mark = ","), length(traits)))

  present_ep <- Filter(function(e) file.exists(src_path(e)), endpoints)
  all_rows <- list(); pairs <- list(); het <- list()
  for (x in Filter(function(t) file.exists(src_path(t)), c(SHARED, CONTROL))) {
    log_path <- run_rg(x, present_ep)
    pr <- parse_run(x, present_ep, log_path)
    all_rows[[length(all_rows) + 1L]] <- pr$rows
    rg <- setNames(as.list(pr$rows$rg), pr$rows$trait2)
    if (x %in% SHARED) for (fam in names(DEFINITION_SETS)) {
      pairs <- c(pairs, pairwise(x, fam, DEFINITION_SETS[[fam]], rg, pr$pseudo))
      h <- heterogeneity(x, fam, DEFINITION_SETS[[fam]], rg, pr$pseudo)
      if (!is.null(h)) het[[length(het) + 1L]] <- h
    }
  }

  ## positive control: MDD_EHR vs MDD_Clin against each other trait (not MDD, which
  ## contains both)
  if (all(vapply(CONTROL, function(t) file.exists(src_path(t)), logical(1)))) {
    others <- c(Filter(function(t) file.exists(src_path(t)), c("ASD", "SCZ", "BIP", "PTSD")), present_ep)
    for (y in others) {
      rg_y <- list(); pv_y <- list()
      if (y %in% present_ep) {
        for (d in CONTROL) {
          pr <- parse_run(d, present_ep, file.path(RUNS, sprintf("rg_%s.log", d)))
          rec <- pr$rows[trait2 == y]
          rg_y[[d]] <- rec$rg; pv_y[[d]] <- pr$pseudo[[y]]
        }
      } else {
        log_path <- run_rg(y, CONTROL, tag = "__control")
        pr <- parse_run(y, CONTROL, log_path)
        all_rows[[length(all_rows) + 1L]] <- pr$rows
        rg_y <- setNames(as.list(pr$rows$rg), pr$rows$trait2)
        pv_y <- pr$pseudo
      }
      if (all(vapply(CONTROL, function(d) !is.null(pv_y[[d]]), logical(1))))
        pairs <- c(pairs, pairwise(y, "MDD definition (positive control)", CONTROL, rg_y, pv_y))
    }
    log_path <- run_rg("MDD_EHR", "MDD_Clin", tag = "__control")
    all_rows[[length(all_rows) + 1L]] <- parse_run("MDD_EHR", "MDD_Clin", log_path)$rows
  }

  rgc <- unique(rbindlist(all_rows, fill = TRUE), by = c("trait1", "trait2"))
  pw  <- if (length(pairs)) rbindlist(pairs) else data.table()
  hz  <- if (length(het))   rbindlist(het)   else data.table()

  dist <- fread(file.path(ROOT, "results", "definitions", "pairwise_distance.tsv"), sep = "\t")
  rule <- as.integer(dist$mode_differs | dist$conditions_differ | dist$controls_differ)
  dist[, distance := 1 - (icd10_jaccard + source_jaccard + (1 - rule)) / 3]
  keep <- c("def1", "def2", "distance", "icd10_jaccard", "source_jaccard", "log_case_ratio")
  rev  <- copy(dist); setnames(rev, c("def1", "def2"), c("def2", "def1"))
  both <- unique(rbind(dist[, ..keep], rev[, ..keep]), by = c("def1", "def2"))
  if (nrow(pw)) pw <- merge(pw, both, by = c("def1", "def2"), all.x = TRUE)

  fwrite(rgc, file.path(OUT, "rg_common.tsv"), sep = "\t")
  fwrite(pw,  file.path(OUT, "pairwise_delta.tsv"), sep = "\t")
  fwrite(hz,  file.path(OUT, "heterogeneity.tsv"), sep = "\t")
  bad <- rgc[!is.na(se_check_ok) & se_check_ok == FALSE]
  if (nrow(bad)) stop(sprintf("%d pairs failed the SE rebuild check", nrow(bad)), call. = FALSE)
  print(hz)
  print(pw)
}

if (sys.nframe() == 0L) main()
