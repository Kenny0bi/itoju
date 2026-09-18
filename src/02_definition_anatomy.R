## Turn FinnGen R12 endpoint definitions into comparable code sets and distances.
##
## R implementation of src/02_definition_anatomy.py. Deterministic throughout, so the
## three output tables match the Python exactly. The .xlsx endpoint file is read with
## readxl rather than pandas; the parsed contents are the same.
##
## For every endpoint in the study, this script:
##   1. resolves INCLUDE chains, so an endpoint like SLEEP inherits the rules of the
##      endpoints it unions;
##   2. expands each registry's ICD-10 patterns against the WHO ICD-10 code universe
##      (the codes present in the phecodeX WHO map, plus their three-character parents);
##   3. records which data sources can make someone a case: hospital discharge (HD),
##      cause of death (COD), primary care (OUTPAT, Avohilmo), drug purchases (KELA_ATC),
##      drug reimbursement (KELA_REIMB), and extra case conditions (CONDITIONS) or
##      control rules;
##   4. maps the expanded ICD-10 set to phecodeX, the shared standard vocabulary;
##   5. computes pairwise definition distances within each condition family.
##
## FinnGen pattern syntax (FinnGen endpoint file format description):
##   "|" separates alternatives, each a regular expression matched from the start of
##   the code; "%" marks a "mode" rule (the code must be the most common among the
##   related diagnoses); "$!$" means the registry is deliberately not used;
##   "!NAME" in CONDITIONS means "and not a case of NAME".
##
## Outputs: results/definitions/endpoint_anatomy.tsv, code_membership.tsv,
## pairwise_distance.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1L])) else getwd()
  normalizePath(file.path(d, ".."), mustWork = FALSE)
}
ROOT        <- find_root()
ENDPOINTS   <- file.path(ROOT, "data", "raw", "endpoints", "finngen_R12_endpoint_core_noncore_1.0.xlsx")
MANIFEST    <- file.path(ROOT, "data", "raw", "finngen_R12_manifest.tsv")
PHECODE_MAP <- file.path(ROOT, "data", "ref", "phecode", "phecodeX_ICD_WHO_map_flat.csv")
OUT         <- file.path(ROOT, "results", "definitions")

FAMILIES <- list(
  Epilepsy = c("G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE", "G6_STATUSEPI"),
  Sleep    = c("SLEEP", "G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "F5_INSOMNIA",
               "KRA_PSY_SLEEP_NONORG_EXMORE", "G6_SLEEPDISOTH", "F5_SLEEP_NOS"),
  Bowel    = c("K11_CONSTIPATION", "K11_OTHFUNC", "K11_IBS", "K11_FUNCDYSP", "K11_REFLUX"),
  Bladder  = c("N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD"),
  ADHD     = c("F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"),
  `Intellectual disability` = c("F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE"),
  Autism   = c("KRA_PSY_DEVWIDE_EXMORE", "KRA_PSY_AUTISM_EXMORE")
)

ICD10_COLS    <- c(HD = "HD_ICD_10", COD = "COD_ICD_10", OUTPAT = "OUTPAT_ICD")
OTHER_SOURCES <- c(KELA_ATC = "KELA_ATC", KELA_REIMB = "KELA_REIMB", OPER = "OPER_NOM")

present <- function(v) !is.na(v) && !(trimws(as.character(v)) %in% c("", "nan", "$!$"))

icd10_universe <- function() {
  m <- fread(PHECODE_MAP, colClasses = "character", encoding = "Latin-1")
  codes <- gsub(".", "", m$icd, fixed = TRUE)
  list(map = m, universe = sort(unique(c(codes, substr(codes, 1L, 3L)))))
}

split_alternatives <- function(pattern) {
  ## Split on "|" outside square brackets. FinnGen writes F84[0|5], where the bar
  ## sits inside a character class, so a plain split would cut the class in half.
  chars <- strsplit(pattern, "", fixed = TRUE)[[1L]]
  parts <- character(0); depth <- 0L; cur <- ""
  for (ch in chars) {
    if (ch == "[") depth <- depth + 1L else if (ch == "]") depth <- depth - 1L
    if (ch == "|" && depth == 0L) { parts <- c(parts, cur); cur <- "" } else cur <- paste0(cur, ch)
  }
  c(parts, cur)
}

expand_pattern <- function(pattern, universe) {
  if (!present(pattern)) return(list(codes = character(0), mode = FALSE))
  out <- character(0); mode <- FALSE
  for (alt in split_alternatives(as.character(pattern))) {
    alt <- trimws(alt)
    if (startsWith(alt, "%")) { mode <- TRUE; alt <- substring(alt, 2L) }
    alt <- gsub(".", "", alt, fixed = TRUE)
    if (!nzchar(alt)) next
    out <- c(out, universe[grepl(paste0("^", alt), universe)])
  }
  list(codes = unique(out), mode = mode)
}

resolve <- function(name, table, universe, seen = character(0)) {
  ## Case rules for an endpoint, with INCLUDE chains folded in.
  if (name %in% seen || !(name %in% rownames(table))) return(NULL)
  seen <- c(seen, name)
  r <- table[name, ]
  rule <- list(icd10 = setNames(vector("list", length(ICD10_COLS)), names(ICD10_COLS)),
               sources = character(0), mode = FALSE, atc = character(0),
               conditions = character(0), includes = character(0))
  for (nm in names(rule$icd10)) rule$icd10[[nm]] <- character(0)
  for (src in names(ICD10_COLS)) {
    e <- expand_pattern(r[[ICD10_COLS[[src]]]], universe)
    if (length(e$codes)) {
      rule$icd10[[src]] <- union(rule$icd10[[src]], e$codes)
      rule$sources <- union(rule$sources, src)
    }
    rule$mode <- rule$mode || e$mode
  }
  for (src in names(OTHER_SOURCES)) {
    v <- r[[OTHER_SOURCES[[src]]]]
    if (present(v)) {
      rule$sources <- union(rule$sources, src)
      if (src == "KELA_ATC") rule$atc <- union(rule$atc, strsplit(as.character(v), "|", fixed = TRUE)[[1L]])
    }
  }
  if (present(r[["CONDITIONS"]])) rule$conditions <- c(rule$conditions, as.character(r[["CONDITIONS"]]))
  if (present(r[["INCLUDE"]])) {
    for (child in trimws(strsplit(as.character(r[["INCLUDE"]]), "|", fixed = TRUE)[[1L]])) {
      rule$includes <- c(rule$includes, child)
      sub <- resolve(child, table, universe, seen)
      if (is.null(sub)) next
      for (s in names(ICD10_COLS)) rule$icd10[[s]] <- union(rule$icd10[[s]], sub$icd10[[s]])
      rule$sources    <- union(rule$sources, sub$sources)
      rule$mode       <- rule$mode || sub$mode
      rule$atc        <- union(rule$atc, sub$atc)
      rule$conditions <- c(rule$conditions, sub$conditions)
    }
  }
  rule
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) NA_real_ else length(intersect(a, b)) / length(u)
}

main <- function() {
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  tbl <- as.data.frame(read_excel(ENDPOINTS))
  rownames(tbl) <- tbl$NAME
  man <- fread(MANIFEST, sep = "\t")
  setkey(man, phenocode)

  pu <- icd10_universe()
  phe_map <- pu$map; universe <- pu$universe
  phe_map[, icd_flat := gsub(".", "", icd, fixed = TRUE)]
  icd_to_phe <- split(phe_map$phecode, phe_map$icd_flat)
  phe_name <- setNames(phe_map$phecode_string, phe_map$phecode)

  rows <- list(); members <- list(); rules <- list()
  for (fam in names(FAMILIES)) for (n in FAMILIES[[fam]]) {
    rule <- resolve(n, tbl, universe)
    rules[[n]] <- rule
    all_icd <- unique(unlist(rule$icd10))
    phecodes <- if (length(all_icd)) unique(unlist(icd_to_phe[all_icd])) else character(0)
    phecodes <- phecodes[!is.na(phecodes)]
    mapped <- all_icd[all_icd %in% names(icd_to_phe)]
    r <- tbl[n, ]
    rows[[length(rows) + 1L]] <- data.table(
      family = fam, endpoint = n, longname = r[["LONGNAME"]],
      cases = as.integer(man[n, num_cases]), controls = as.integer(man[n, num_controls]),
      n_icd10 = length(all_icd), icd10_codes = paste(sort(all_icd), collapse = " "),
      sources = paste(sort(rule$sources), collapse = " "), mode_rule = rule$mode,
      atc = paste(sort(rule$atc), collapse = " "), includes = paste(rule$includes, collapse = " "),
      case_conditions = paste(rule$conditions, collapse = " ; "),
      control_exclude = if (present(r[["CONTROL_EXCLUDE"]])) as.character(r[["CONTROL_EXCLUDE"]]) else "",
      control_conditions = if (present(r[["CONTROL_CONDITIONS"]])) as.character(r[["CONTROL_CONDITIONS"]]) else "",
      icd9_used = present(r[["HD_ICD_9"]]), icd8_used = present(r[["HD_ICD_8"]]),
      phecodes = paste(sort(phecodes), collapse = " "),
      phecode_names = paste(sort(unname(phe_name[phecodes])), collapse = " ; "),
      icd10_mapped_frac = if (length(all_icd)) length(mapped) / length(all_icd) else NA_real_)
    for (src in names(rule$icd10)) for (c in rule$icd10[[src]])
      members[[length(members) + 1L]] <- data.table(family = fam, endpoint = n, source = src, icd10 = c)
  }
  anatomy <- rbindlist(rows)
  fwrite(anatomy, file.path(OUT, "endpoint_anatomy.tsv"), sep = "\t")
  fwrite(rbindlist(members), file.path(OUT, "code_membership.tsv"), sep = "\t")

  setkey(anatomy, endpoint)
  pairs <- list()
  for (fam in names(FAMILIES)) {
    names_f <- FAMILIES[[fam]]
    if (length(names_f) < 2L) next
    cb <- combn(names_f, 2L)
    for (i in seq_len(ncol(cb))) {
      x <- cb[1L, i]; yy <- cb[2L, i]
      ax <- anatomy[x]; ay <- anatomy[yy]
      cx <- strsplit(ax$icd10_codes, " ")[[1L]]; cy <- strsplit(ay$icd10_codes, " ")[[1L]]
      sx <- strsplit(ax$sources, " ")[[1L]];     sy <- strsplit(ay$sources, " ")[[1L]]
      px <- strsplit(ax$phecodes, " ")[[1L]];    py <- strsplit(ay$phecodes, " ")[[1L]]
      ctrl_same <- ax$control_exclude == ay$control_exclude &&
                   ax$control_conditions == ay$control_conditions
      pairs[[length(pairs) + 1L]] <- data.table(
        family = fam, def1 = x, def2 = yy,
        icd10_jaccard = jaccard(cx, cy), source_jaccard = jaccard(sx, sy),
        phecode_jaccard = jaccard(px, py),
        mode_differs = ax$mode_rule != ay$mode_rule,
        conditions_differ = ax$case_conditions != ay$case_conditions,
        controls_differ = !ctrl_same,
        log_case_ratio = abs(log(ax$cases / ay$cases)))
    }
  }
  pw <- rbindlist(pairs)
  fwrite(pw, file.path(OUT, "pairwise_distance.tsv"), sep = "\t")

  print(anatomy[, .(family, endpoint, cases, n_icd10, sources, mode_rule, atc,
                    control_exclude, icd10_mapped_frac, phecode_names)])
  cat("\n")
  print(pw)
}

if (sys.nframe() == 0L) main()
