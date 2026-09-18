## Stream one FinnGen R12 endpoint, keep HapMap3 SNPs, and write an LDSC-ready file.
##
## R implementation of src/01_fetch_munge_finngen.py. Same filtering rule, same
## output columns, same resumption markers, and the same md5 of the streamed gzip
## recorded for provenance, so either language can produce the munged inputs.
##
## The raw FinnGen file (0.3 to 0.8 GB) is never written to disk. It is read straight
## from the public bucket, filtered to the HapMap3 rsIDs in the LDSC w_hm3.snplist,
## and saved as a small tab-separated file that munge_sumstats.py then turns into
## .sumstats.gz.
##
## FinnGen R12 columns: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta
## af_alt af_alt_cases af_alt_controls (GRCh38; beta is for alt). LDSC merges on rsID,
## so the genome build does not matter for the regression itself.
##
## munge_sumstats.py is Python and has no R equivalent; this script shells out to it
## exactly as the Python version does, so both languages call the same tool.
##
## Usage: Rscript src/01_fetch_munge_finngen.R ENDPOINT [ENDPOINT ...]
## Resumable: an endpoint with a finished marker in data/munged is skipped.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT     <- find_root()
TOOLS    <- normalizePath(file.path(ROOT, "..", "tools"), mustWork = FALSE)
PY       <- file.path(TOOLS, "ldsc_env", "bin", "python")
MUNGE    <- file.path(TOOLS, "cbiit_ldsc", "munge_sumstats.py")
HM3      <- file.path(ROOT, "data", "ref", "z8182036", "eur_w_ld_chr", "w_hm3.snplist")
MANIFEST <- file.path(ROOT, "data", "raw", "finngen_R12_manifest.tsv")
OUT      <- file.path(ROOT, "data", "munged")
TMP      <- file.path(ROOT, "data", "processed", "hm3_filtered")

MAF_MIN <- 0.01

load_manifest <- function() {
  m <- fread(MANIFEST, sep = "\t", colClasses = "character")
  split(m, m$phenocode)
}

load_hm3 <- function() {
  x <- fread(HM3, sep = "\t", select = 1L, colClasses = "character")
  unique(x[[1L]])
}

stream_filter <- function(url, hm3, dest, keys_path = NULL) {
  ## Download and filter in one pass. Returns rows read, rows kept, and the md5 of
  ## the gzip stream, so the exact bytes processed are recorded.
  local_gz <- tempfile(fileext = ".gz")
  on.exit(unlink(local_gz), add = TRUE)
  utils::download.file(url, local_gz, mode = "wb", quiet = TRUE)
  md5 <- digest(file = local_gz, algo = "md5")

  con <- gzfile(local_gz, "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(sub("^#", "", readLines(con, n = 1L)), "\t", fixed = TRUE)[[1L]]
  ix <- setNames(seq_along(header), header)

  out_con <- file(dest, "wt"); on.exit(close(out_con), add = TRUE)
  writeLines("SNP\tA1\tA2\tBETA\tSE\tP\tFRQ", out_con)
  keys_con <- if (!is.null(keys_path)) file(paste0(keys_path, ".partial"), "wt") else NULL
  if (!is.null(keys_con)) on.exit(close(keys_con), add = TRUE)

  n_read <- 0L; n_kept <- 0L
  repeat {
    chunk <- readLines(con, n = 200000L)
    if (!length(chunk)) break
    n_read <- n_read + length(chunk)
    f <- tstrsplit(chunk, "\t", fixed = TRUE)
    rsids <- f[[ix[["rsids"]]]]
    keep  <- nzchar(rsids)
    if (!any(keep)) next
    ## a position can carry several rsIDs separated by commas
    hit <- vapply(strsplit(rsids, ",", fixed = TRUE),
                  function(v) { m <- v[v %in% hm3]; if (length(m)) m[1L] else NA_character_ },
                  character(1L))
    af  <- suppressWarnings(as.numeric(f[[ix[["af_alt"]]]]))
    keep <- keep & !is.na(hit) & !is.na(af) & pmin(af, 1 - af) >= MAF_MIN
    if (!any(keep)) next
    writeLines(paste(hit[keep], f[[ix[["alt"]]]][keep], f[[ix[["ref"]]]][keep],
                     f[[ix[["beta"]]]][keep], f[[ix[["sebeta"]]]][keep],
                     f[[ix[["pval"]]]][keep], af[keep], sep = "\t"), out_con)
    if (!is.null(keys_con))
      writeLines(paste(paste(f[[ix[["chrom"]]]][keep], f[[ix[["pos"]]]][keep],
                             f[[ix[["ref"]]]][keep], f[[ix[["alt"]]]][keep], sep = ":"),
                       hit[keep], sep = "\t"), keys_con)
    n_kept <- n_kept + sum(keep)
  }
  list(n_read = n_read, n_kept = n_kept, md5 = md5)
}

run_endpoint <- function(endpoint, manifest, hm3) {
  done <- file.path(OUT, paste0(endpoint, ".done"))
  if (file.exists(done)) { cat(endpoint, ": already done\n", sep = ""); return(invisible(NULL)) }
  row <- manifest[[endpoint]]
  dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(TMP, paste0(endpoint, ".hm3.tsv"))
  t0  <- Sys.time()

  ## Once, save the FinnGen GRCh38 variant key -> rsID map for HapMap3 SNPs. The
  ## FinnGen LD matrix names variants as chr:pos:ref:alt, so this map is what turns
  ## it into rsID-keyed LD scores.
  keys_path <- file.path(ROOT, "data", "processed", "finngen_hm3_variant_keys.tsv")
  want_keys <- if (file.exists(keys_path)) NULL else keys_path

  r <- stream_filter(row$path_https, hm3, tmp, want_keys)
  if (!is.null(want_keys)) file.rename(paste0(keys_path, ".partial"), keys_path)

  cat(sprintf("%s: read %s variants, kept %s HapMap3 with MAF >= %g, md5 %s, %.0fs\n",
              endpoint, format(r$n_read, big.mark = ","), format(r$n_kept, big.mark = ","),
              MAF_MIN, r$md5, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  system2(PY, c(MUNGE, "--sumstats", tmp,
                "--N-cas", row$num_cases, "--N-con", row$num_controls,
                "--signed-sumstats", "BETA,0", "--merge-alleles", HM3,
                "--out", file.path(OUT, endpoint)), stdout = FALSE)
  unlink(tmp)
  writeLines(c(paste0("endpoint\t", endpoint), paste0("cases\t", row$num_cases),
               paste0("controls\t", row$num_controls), paste0("url\t", row$path_https),
               paste0("md5_gz\t", r$md5), paste0("variants_read\t", r$n_read),
               paste0("hm3_kept\t", r$n_kept),
               paste0("finished\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), done)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  manifest <- load_manifest()
  hm3 <- load_hm3()
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  for (e in args) run_endpoint(e, manifest, hm3)
}

if (sys.nframe() == 0L) main()
