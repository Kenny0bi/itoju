## SNP heritability of every munged trait, and genetic correlations for the design.
## rg runs carry no prevalence flags: rg is the same on either scale.
##
## Usage: Rscript src/03_ldsc.R

suppressPackageStartupMessages(library(data.table))

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
OUT    <- file.path(ROOT, "results", "ldsc")

## trait: c(sample prevalence, population prevalence); NULL = report observed scale only
PSYCH <- list(
  ASD      = c(18381 / (18381 + 27969), 0.012),  # Grove 2019, N cases and controls, K from the paper
  SCZ      = c(0.5, 0.01),
  BIP      = c(0.5, 0.02),
  MDD      = c(0.5, 0.15),
  PTSD     = NULL,                               # z-score meta-analysis mixing symptom scores and diagnoses
  MDD_EHR  = c(0.5, 0.15),
  MDD_Clin = c(0.5, 0.15)
)
SUMSTATS <- c(ASD = "ASD_Grove2019")

ENV <- c("OMP_NUM_THREADS=2", "OPENBLAS_NUM_THREADS=2", "MKL_NUM_THREADS=2", "VECLIB_MAXIMUM_THREADS=2")
NUM <- "(-?[0-9.]+(?:e-?[0-9]+)?)"

trait_path <- function(trait) {
  nm <- if (trait %in% names(SUMSTATS)) SUMSTATS[[trait]] else trait
  file.path(MUNGED, paste0(nm, ".sumstats.gz"))
}

ldsc_run <- function(args, out, done_marker) {
  log_path <- paste0(out, ".log")
  if (file.exists(log_path)) {
    txt <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
    if (grepl(done_marker, txt, fixed = TRUE)) return(txt)
  }
  system2("nice", c("-n", "10", PY, LDSC, args, "--ref-ld-chr", LD, "--w-ld-chr", LD, "--out", out),
          stdout = FALSE, stderr = FALSE, env = ENV)
  txt <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  if (!grepl(done_marker, txt, fixed = TRUE))
    stop(sprintf("LDSC produced no result for %s; see %s", basename(out), log_path), call. = FALSE)
  txt
}

grab <- function(text, label) {
  m <- regmatches(text, regexpr(paste0(gsub("([\\^$.|?*+()\\[\\]{}])", "\\\\\\1", label, perl = TRUE),
                                       ":\\s*", NUM, "(?:\\s*\\(", NUM, "\\))?"), text, perl = TRUE))
  if (!length(m)) return(c(NA_real_, NA_real_))
  nums <- as.numeric(regmatches(m, gregexpr(NUM, m, perl = TRUE))[[1L]])
  c(if (length(nums) >= 1L) nums[1L] else NA_real_,
    if (length(nums) >= 2L) nums[2L] else NA_real_)
}

h2_record <- function(trait, text, kind, prev) {
  scale <- if (!is.null(prev)) "Liability" else "Observed"
  h2    <- grab(text, sprintf("Total %s scale h2", scale))
  icpt  <- grab(text, "Intercept")
  data.table(trait = trait, kind = kind,
             samp_prev = if (!is.null(prev)) prev[1L] else NA_real_,
             pop_prev  = if (!is.null(prev)) prev[2L] else NA_real_,
             scale = tolower(scale), h2 = h2[1L], h2_se = h2[2L],
             h2_z = if (!is.na(h2[1L]) && !is.na(h2[2L]) && h2[2L] != 0) h2[1L] / h2[2L] else NA_real_,
             lambda_gc = grab(text, "Lambda GC")[1L], mean_chi2 = grab(text, "Mean Chi^2")[1L],
             intercept = icpt[1L], intercept_se = icpt[2L])
}

rg_record <- function(t1, t2, text) {
  tail_lines <- strsplit(substring(text, regexpr("Summary of Genetic Correlation Results", text)), "\n")[[1L]]
  hdr <- strsplit(trimws(tail_lines[2L]), "\\s+")[[1L]]
  val <- strsplit(trimws(tail_lines[3L]), "\\s+")[[1L]]
  row <- list(trait1 = t1, trait2 = t2)
  for (i in seq_along(hdr)) {
    if (hdr[i] %in% c("p1", "p2")) next
    row[[hdr[i]]] <- if (val[i] == "NA") NA_real_ else as.numeric(val[i])
  }
  as.data.table(row)
}

main <- function() {
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  psych <- names(PSYCH)[vapply(names(PSYCH), function(t) file.exists(trait_path(t)), logical(1))]

  endpoints <- list()
  for (done in sort(list.files(MUNGED, pattern = "\\.done$", full.names = TRUE))) {
    kv <- do.call(rbind, strsplit(readLines(done, warn = FALSE), "\t", fixed = TRUE))
    meta <- setNames(as.list(kv[, 2L]), kv[, 1L])
    if (!is.null(meta$endpoint)) {
      ca <- as.numeric(meta$cases); co <- as.numeric(meta$controls)
      endpoints[[length(endpoints) + 1L]] <- list(name = meta$endpoint, prev = ca / (ca + co))
    }
  }

  h2_rows <- list()
  for (t in psych) {
    prev <- PSYCH[[t]]
    args <- c("--h2", trait_path(t))
    if (!is.null(prev)) args <- c(args, "--samp-prev", prev[1L], "--pop-prev", prev[2L])
    nm   <- if (t %in% names(SUMSTATS)) SUMSTATS[[t]] else t
    text <- ldsc_run(args, file.path(OUT, paste0("h2_", nm)), "scale h2")
    h2_rows[[length(h2_rows) + 1L]] <- h2_record(t, text, "psychiatric", prev)
  }
  for (e in endpoints) {
    text <- ldsc_run(c("--h2", trait_path(e$name), "--samp-prev", e$prev, "--pop-prev", e$prev),
                     file.path(OUT, paste0("h2_", e$name)), "Total Liability scale h2")
    h2_rows[[length(h2_rows) + 1L]] <- h2_record(e$name, text, "finngen", c(e$prev, e$prev))
  }

  rg_rows <- list()
  for (t in psych) for (e in endpoints) {
    text <- ldsc_run(c("--rg", paste(trait_path(t), trait_path(e$name), sep = ",")),
                     file.path(OUT, sprintf("rg_%s_%s", t, e$name)),
                     "Summary of Genetic Correlation Results")
    rg_rows[[length(rg_rows) + 1L]] <- rg_record(t, e$name, text)
  }
  if (length(psych) >= 2L) {
    cb <- combn(psych, 2L)
    for (i in seq_len(ncol(cb))) {
      t1 <- cb[1L, i]; t2 <- cb[2L, i]
      text <- ldsc_run(c("--rg", paste(trait_path(t1), trait_path(t2), sep = ",")),
                       file.path(OUT, sprintf("rg_%s_%s", t1, t2)),
                       "Summary of Genetic Correlation Results")
      rg_rows[[length(rg_rows) + 1L]] <- rg_record(t1, t2, text)
    }
  }

  h2 <- rbindlist(h2_rows, fill = TRUE)
  rg <- rbindlist(rg_rows, fill = TRUE)
  fwrite(h2, file.path(OUT, "h2.tsv"), sep = "\t")
  fwrite(rg, file.path(OUT, "rg.tsv"), sep = "\t")
  print(h2)
  print(rg)
}

if (sys.nframe() == 0L) main()
