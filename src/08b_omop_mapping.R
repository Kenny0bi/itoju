## Map every FinnGen endpoint definition to the OMOP standard vocabularies.
##
## R implementation of src/08b_omop_mapping.py. Deterministic: same pattern expansion,
## same "Maps to" traversal, same ancestor closure depth, so the three output tables
## match the Python.
##
## Why: FinnGen writes definitions as regular expressions over Finnish ICD-10, ICD-9
## and ICD-8 codes plus Finnish drug and reimbursement codes. To compare definitions
## with each other, and to rebuild them in any OMOP CDM database, each definition is
## translated to standard concepts:
##   - ICD-10 (WHO) codes -> "Maps to" -> SNOMED CT standard condition concepts;
##   - Finnish ICD-9 codes -> ICD-9-CM (approximate: Finnish ICD-9 uses 4 digits plus
##     an optional letter, e.g. 5965B; the letter is dropped and the code is dotted,
##     596.5) -> SNOMED;
##   - ICD-8 has no OMOP vocabulary, so ICD-8 rules are counted as unmappable;
##   - ATC classes (e.g. A06A, drugs for constipation) -> RxNorm ingredients via
##     CONCEPT_ANCESTOR.
##
## Overlap between two definitions of the same condition is then measured three ways,
## from coarse to clinical: ICD-10 code Jaccard (from 02), phecodeX Jaccard (from 02),
## SNOMED exact-concept Jaccard, and SNOMED hierarchy-aware Jaccard (each concept
## expanded to its ancestors up to MAX_LEVELS levels, so a parent and a child concept
## count as related rather than disjoint).
##
## Inputs: data/ref/omop/* (src/08a_extract_omop.sh), results/definitions/
## endpoint_anatomy.tsv and pairwise_distance.tsv (src/02), the R12 endpoint file.
## Outputs: results/omop/endpoint_omop.tsv, endpoint_snomed_concepts.tsv,
## pairwise_omop.tsv

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
ROOT      <- find_root()
OMOP      <- file.path(ROOT, "data", "ref", "omop")
DEF       <- file.path(ROOT, "results", "definitions")
ENDPOINTS <- file.path(ROOT, "data", "raw", "endpoints", "finngen_R12_endpoint_core_noncore_1.0.xlsx")
OUT       <- file.path(ROOT, "results", "omop")
MAX_LEVELS <- 2L   # hierarchy-aware overlap: include ancestors within two levels

ICD9_COLS <- c("HD_ICD_9", "COD_ICD_9")
ICD8_COLS <- c("HD_ICD_8", "COD_ICD_8")

present <- function(v) !is.na(v) && !(trimws(as.character(v)) %in% c("", "nan", "$!$"))

split_alternatives <- function(pattern) {
  chars <- strsplit(pattern, "", fixed = TRUE)[[1L]]
  parts <- character(0); depth <- 0L; cur <- ""
  for (ch in chars) {
    if (ch == "[") depth <- depth + 1L else if (ch == "]") depth <- depth - 1L
    if (ch == "|" && depth == 0L) { parts <- c(parts, cur); cur <- "" } else cur <- paste0(cur, ch)
  }
  c(parts, cur)
}

dot_icd10 <- function(code) ifelse(nchar(code) <= 3L, code,
                                   paste0(substr(code, 1L, 3L), ".", substring(code, 4L)))

read_omop <- function(name) fread(file.path(OMOP, name), sep = "\t", colClasses = "character",
                                  quote = "", na.strings = NULL)

icd9_patterns <- function(name, table, seen = character(0)) {
  ## Finnish ICD-9 patterns for an endpoint, following INCLUDE chains; plus whether
  ## ICD-8 is used.
  if (name %in% seen || !(name %in% rownames(table))) return(list(pats = character(0), icd8 = FALSE))
  seen <- c(seen, name)
  r <- table[name, ]
  pats <- character(0)
  for (c in ICD9_COLS) if (present(r[[c]]))
    pats <- c(pats, sub("^%", "", trimws(split_alternatives(as.character(r[[c]])))))
  icd8 <- any(vapply(ICD8_COLS, function(c) present(r[[c]]), logical(1)))
  if (present(r[["INCLUDE"]])) {
    for (child in trimws(strsplit(as.character(r[["INCLUDE"]]), "|", fixed = TRUE)[[1L]])) {
      sub <- icd9_patterns(child, table, seen)
      pats <- c(pats, sub$pats); icd8 <- icd8 || sub$icd8
    }
  }
  list(pats = pats, icd8 = icd8)
}

main <- function() {
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  src   <- read_omop("source_concepts.tsv")
  maps  <- read_omop("maps_to.tsv")
  tgt   <- read_omop("target_concepts.tsv"); setkey(tgt, concept_id)
  anc   <- read_omop("target_ancestors.tsv")
  atc_desc <- read_omop("atc_descendants.tsv")
  related  <- read_omop("related_concepts.tsv"); setkey(related, concept_id)
  versions <- if (file.exists(file.path(OMOP, "vocabulary_versions.tsv"))) read_omop("vocabulary_versions.tsv") else NULL

  maps_to <- split(maps$concept_id_2, maps$concept_id_1)
  icd10   <- setNames(src[vocabulary_id == "ICD10", concept_id], src[vocabulary_id == "ICD10", concept_code])
  icd9    <- src[vocabulary_id == "ICD9CM"]
  icd9_undotted <- setNames(icd9$concept_id, gsub(".", "", icd9$concept_code, fixed = TRUE))
  atc     <- setNames(src[vocabulary_id == "ATC", concept_id], src[vocabulary_id == "ATC", concept_code])

  anc <- anc[as.integer(min_levels_of_separation) <= MAX_LEVELS]
  closure <- split(anc$ancestor_concept_id, anc$descendant_concept_id)
  snomed_std <- tgt[vocabulary_id == "SNOMED" & standard_concept == "S", concept_id]

  anatomy <- fread(file.path(DEF, "endpoint_anatomy.tsv"), sep = "\t")
  tbl <- as.data.frame(read_excel(ENDPOINTS)); rownames(tbl) <- tbl$NAME

  rows <- list(); concept_rows <- list(); sets <- list()
  for (i in seq_len(nrow(anatomy))) {
    a <- anatomy[i]
    e <- a$endpoint
    codes <- if (!is.na(a$icd10_codes) && nzchar(a$icd10_codes)) strsplit(a$icd10_codes, " ")[[1L]] else character(0)
    cids  <- unname(icd10[dot_icd10(codes)])
    found <- cids[!is.na(cids)]
    snomed <- if (length(found)) intersect(unique(unlist(maps_to[found])), snomed_std) else character(0)

    ip <- icd9_patterns(e, tbl)
    icd9_hits <- character(0)
    for (p in ip$pats) {
      rx <- paste0("^", sub("[A-Z]$", "", gsub(".", "", p, fixed = TRUE)))
      icd9_hits <- union(icd9_hits, unname(icd9_undotted[grepl(rx, names(icd9_undotted))]))
    }
    snomed9 <- if (length(icd9_hits)) intersect(unique(unlist(maps_to[icd9_hits])), snomed_std) else character(0)

    atc_codes <- if (!is.na(a$atc) && nzchar(trimws(a$atc))) strsplit(a$atc, " ")[[1L]] else character(0)
    ingredients <- character(0)
    for (code in atc_codes) {
      cid <- atc[[code]]
      if (is.null(cid) || is.na(cid)) next
      desc <- atc_desc[ancestor_concept_id == cid, descendant_concept_id]
      if (!length(desc)) next
      hit <- related[J(desc), nomatch = 0L][vocabulary_id == "RxNorm" & concept_class_id == "Ingredient", concept_id]
      ingredients <- union(ingredients, hit)
    }

    all_snomed <- union(snomed, snomed9)
    sets[[e]] <- all_snomed
    rows[[length(rows) + 1L]] <- data.table(
      family = a$family, endpoint = e, cases = a$cases,
      icd10_codes = length(codes), icd10_in_omop = length(found),
      icd10_mapped_to_snomed = sum(vapply(found, function(c) length(intersect(maps_to[[c]], snomed_std)) > 0L, logical(1))),
      snomed_from_icd10 = length(snomed),
      icd9_patterns = length(ip$pats), icd9cm_concepts_matched = length(icd9_hits),
      snomed_from_icd9 = length(snomed9),
      snomed_total = length(all_snomed), uses_icd8_unmappable = ip$icd8,
      atc_classes = paste(atc_codes, collapse = " "), rxnorm_ingredients = length(ingredients))
    for (c in sort(all_snomed))
      concept_rows[[length(concept_rows) + 1L]] <- data.table(
        endpoint = e, snomed_concept_id = c, from_icd10 = c %in% snomed,
        domain_id = tgt[J(c), domain_id])
  }

  omop <- rbindlist(rows)
  fwrite(omop, file.path(OUT, "endpoint_omop.tsv"), sep = "\t")
  fwrite(rbindlist(concept_rows), file.path(OUT, "endpoint_snomed_concepts.tsv"), sep = "\t")

  jac <- function(p, q) { u <- union(p, q); if (!length(u)) NA_real_ else length(intersect(p, q)) / length(u) }
  dist <- fread(file.path(DEF, "pairwise_distance.tsv"), sep = "\t")
  out <- list()
  for (i in seq_len(nrow(dist))) {
    d <- dist[i]
    x <- if (!is.null(sets[[d$def1]])) sets[[d$def1]] else character(0)
    y <- if (!is.null(sets[[d$def2]])) sets[[d$def2]] else character(0)
    cx <- union(x, unlist(closure[x])); cy <- union(y, unlist(closure[y]))
    out[[length(out) + 1L]] <- cbind(d, data.table(
      snomed_jaccard = jac(x, y), snomed_hier_jaccard = jac(cx, cy),
      snomed_n1 = length(x), snomed_n2 = length(y)))
  }
  pw <- rbindlist(out)
  fwrite(pw, file.path(OUT, "pairwise_omop.tsv"), sep = "\t")

  if (!is.null(versions)) print(versions)
  print(omop)
  print(pw[, .(family, def1, def2, icd10_jaccard, phecode_jaccard, snomed_jaccard, snomed_hier_jaccard)])
}

if (sys.nframe() == 0L) main()
