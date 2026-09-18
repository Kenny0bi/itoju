## Verify every reference in a BibTeX file against Crossref, and every citation
## against the file. Mismatches fail loudly.
##
## Usage: Rscript src/verify_refs.R paper/references.bib [paper/itoju.tex]

suppressPackageStartupMessages({
  library(jsonlite)
  library(curl)
})

MAILTO <- "itoju-reference-check"

norm_words <- function(s) {
  s <- iconv(if (is.null(s)) "" else s, to = "ASCII//TRANSLIT")
  s <- gsub("\\\\[a-zA-Z]+|[{}\\\\$]", "", s)
  w <- strsplit(trimws(gsub("[^a-z0-9 ]+", " ", tolower(s))), "\\s+")[[1L]]
  w[nzchar(w)]
}

similar <- function(a, b) {
  x <- paste(norm_words(a), collapse = " ")
  y <- paste(norm_words(b), collapse = " ")
  if (!nzchar(x) && !nzchar(y)) return(1)
  1 - adist(x, y)[1, 1] / max(nchar(x), nchar(y), 1L)
}

parse_bib <- function(text) {
  entries <- list()
  ms <- gregexpr("(?s)@(\\w+)\\s*\\{\\s*([^,]+),(.*?)\\n\\}", text, perl = TRUE)[[1L]]
  if (ms[1L] == -1L) return(entries)
  lens <- attr(ms, "match.length")
  for (i in seq_along(ms)) {
    blk  <- substr(text, ms[i], ms[i] + lens[i] - 1L)
    kind <- tolower(sub("(?s)^@(\\w+).*$", "\\1", blk, perl = TRUE))
    key  <- trimws(sub("(?s)^@\\w+\\s*\\{\\s*([^,]+),.*$", "\\1", blk, perl = TRUE))
    fm   <- gregexpr("(\\w+)\\s*=\\s*(\\{(?:[^{}]|\\{(?:[^{}]|\\{[^{}]*\\})*\\})*\\}|\"[^\"]*\"|\\d+)",
                     blk, perl = TRUE)[[1L]]
    flen <- attr(fm, "match.length")
    fields <- list(type = kind)
    if (fm[1L] != -1L) for (j in seq_along(fm)) {
      piece <- substr(blk, fm[j], fm[j] + flen[j] - 1L)
      nm    <- tolower(sub("^(\\w+)\\s*=.*", "\\1", piece, perl = TRUE))
      val   <- trimws(sub("^\\w+\\s*=\\s*", "", piece))
      if (substr(val, 1L, 1L) %in% c("{", "\"")) val <- substr(val, 2L, nchar(val) - 1L)
      fields[[nm]] <- trimws(gsub("\\s+", " ", val))
    }
    entries[[key]] <- fields
  }
  entries
}

fetch_json <- function(url) {
  h <- new_handle(); handle_setheaders(h, "User-Agent" = sprintf("verify_refs/1.0 (%s)", MAILTO))
  r <- tryCatch(curl_fetch_memory(url, handle = h), error = function(e) NULL)
  if (is.null(r)) return(list(status = 0L, body = NULL))
  list(status = r$status_code,
       body = tryCatch(fromJSON(rawToChar(r$content), simplifyVector = FALSE), error = function(e) NULL))
}

crossref <- function(doi) {
  ## Crossref record, or a DataCite record mapped to the same fields (repository
  ## DOIs such as Stanford Digital Repository preprints are registered with
  ## DataCite, not Crossref).
  r <- fetch_json(paste0("https://api.crossref.org/works/", URLencode(doi, reserved = TRUE)))
  if (r$status == 200L && !is.null(r$body)) return(r$body$message)
  if (r$status != 404L && r$status != 0L)
    stop(sprintf("Crossref returned status %s for %s", r$status, doi), call. = FALSE)
  d <- fetch_json(paste0("https://api.datacite.org/dois/", URLencode(doi, reserved = TRUE)))
  if (d$status != 200L || is.null(d$body)) stop(sprintf("neither Crossref nor DataCite resolved %s", doi), call. = FALSE)
  a   <- d$body$data$attributes
  fam <- function(c) if (!is.null(c$familyName)) c$familyName else strsplit(c$name, ",")[[1L]][1L]
  list(title = list(a$titles[[1L]]$title),
       author = lapply(seq_along(a$creators), function(i)
         list(family = fam(a$creators[[i]]), sequence = if (i == 1L) "first" else "additional")),
       issued = list(`date-parts` = list(list(a$publicationYear))),
       `container-title` = list(if (is.null(a$publisher)) "" else a$publisher))
}

first_family <- function(authors_field) {
  first <- trimws(strsplit(authors_field, " and ", fixed = TRUE)[[1L]][1L])
  fam <- if (grepl(",", first, fixed = TRUE)) strsplit(first, ",")[[1L]][1L]
         else tail(strsplit(first, "\\s+")[[1L]], 1L)
  paste(norm_words(fam), collapse = " ")
}

check_entry <- function(key, e) {
  doi <- e$doi
  if (is.null(doi) || !nzchar(doi)) return(c("NO-DOI", "no DOI; check by hand"))
  m <- tryCatch(crossref(doi), error = function(ex) ex)
  if (inherits(m, "condition")) return(c("FAIL", sprintf("Crossref lookup failed: %s", conditionMessage(m))))
  problems <- character(0)
  title <- if (length(m$title)) m$title[[1L]] else ""
  if (similar(if (is.null(e$title)) "" else e$title, title) < 0.90)
    problems <- c(problems, sprintf("title: bib '%s' vs crossref '%s'",
                                    substr(if (is.null(e$title)) "" else e$title, 1, 70), substr(title, 1, 70)))
  years <- unique(unlist(lapply(c("published-print", "published-online", "issued"), function(k)
    if (!is.null(m[[k]]$`date-parts`)) as.character(m[[k]]$`date-parts`[[1L]][[1L]]))))
  if (!is.null(e$year) && nzchar(e$year) && !(e$year %in% years))
    problems <- c(problems, sprintf("year: bib %s vs crossref %s", e$year, paste(sort(years), collapse = ",")))
  authors <- if (is.null(m$author)) list() else m$author
  ## Crossref sometimes lists a consortium (a "name" with no "family") ahead of the
  ## first person, as for Grove et al. 2019. Compare against the first named person,
  ## and also accept a match with the listed first author, so a consortium-first
  ## record neither hides nor fakes a mismatch.
  people <- unlist(lapply(authors, function(a) a$family))
  listed_first <- NULL
  for (a in authors) if (identical(a$sequence, "first")) { listed_first <- if (!is.null(a$family)) a$family else a$name; break }
  if (is.null(listed_first) && length(authors))
    listed_first <- if (!is.null(authors[[1L]]$family)) authors[[1L]]$family else authors[[1L]]$name
  cr_first <- if (length(people)) people[1L] else listed_first
  if (!is.null(e$author) && nzchar(e$author) && (length(cr_first) || length(listed_first))) {
    bf    <- first_family(e$author)
    cands <- unique(vapply(Filter(Negate(is.null), list(cr_first, listed_first)),
                           function(x) paste(norm_words(x), collapse = " "), character(1L)))
    if (!any(vapply(cands, function(c) bf == c || grepl(bf, c, fixed = TRUE) || grepl(c, bf, fixed = TRUE), logical(1))))
      problems <- c(problems, sprintf("first author: bib '%s' vs crossref %s", bf, paste(sort(cands), collapse = ",")))
  }
  container     <- paste(unlist(m$`container-title`), collapse = " ")
  bib_container <- if (!is.null(e$journal)) e$journal else if (!is.null(e$booktitle)) e$booktitle else ""
  if (nzchar(bib_container) && nzchar(container) && similar(bib_container, container) < 0.6 &&
      !all(norm_words(bib_container) %in% norm_words(container)))
    problems <- c(problems, sprintf("journal: bib '%s' vs crossref '%s'", bib_container, container))
  Sys.sleep(0.2)
  if (length(problems)) c("FAIL", paste(problems, collapse = "; "))
  else c("OK", sprintf("%s %s %s", cr_first, paste(sort(years), collapse = ","), substr(container, 1, 40)))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("usage: Rscript src/verify_refs.R references.bib [paper.tex]", call. = FALSE)
  entries  <- parse_bib(paste(readLines(args[1L], warn = FALSE), collapse = "\n"))
  failures <- 0L
  for (key in names(entries)) {
    r <- check_entry(key, entries[[key]])
    failures <- failures + (r[1L] == "FAIL")
    cat(sprintf("%-6s %-28s %s\n", r[1L], key, r[2L]))
  }
  if (length(args) >= 2L) {
    tex    <- paste(readLines(args[2L], warn = FALSE), collapse = "\n")
    groups <- regmatches(tex, gregexpr("\\\\cite[a-z]*\\{[^}]*\\}", tex, perl = TRUE))[[1L]]
    cited  <- unique(trimws(unlist(strsplit(sub("^\\\\cite[a-z]*\\{", "", sub("\\}$", "", groups)), ","))))
    cited  <- cited[nzchar(cited)]
    missing <- sort(setdiff(cited, names(entries)))
    uncited <- sort(setdiff(names(entries), cited))
    for (k in missing) cat(sprintf("FAIL   %-28s cited in tex but not in bib\n", k))
    for (k in uncited) cat(sprintf("FAIL   %-28s in bib but never cited\n", k))
    failures <- failures + length(missing) + length(uncited)
  }
  cat(sprintf("\n%d entries, %d failures\n", length(entries), failures))
  quit(status = if (failures) 1L else 0L)
}

if (sys.nframe() == 0L) main()
