#!/usr/bin/env Rscript

## Shared helpers for the ATCC 24927 GPCR expression analysis.
## The project intentionally keeps raw GEO files untouched and writes derived
## tables under atcc_gpcr_analysis/results.

options(stringsAsFactors = FALSE, warn = 1)

project_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(here) == "atcc_gpcr_analysis") return(here)
  candidate <- file.path(here, "atcc_gpcr_analysis")
  if (dir.exists(candidate)) return(normalizePath(candidate, winslash = "/"))
  normalizePath(file.path(here, ".."), winslash = "/", mustWork = FALSE)
}

PROJ <- project_root()
REPO <- normalizePath(file.path(PROJ, ".."), winslash = "/", mustWork = FALSE)
dir.create(file.path(PROJ, "logs"), recursive = TRUE, showWarnings = FALSE)

start_log <- function(step, inputs = character()) {
  log_path <- file.path(PROJ, "logs", paste0(step, ".log"))
  sink(log_path, split = TRUE)
  cat("ATCC 24927 GPCR analysis\n")
  cat("step:", step, "\n")
  cat("started_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
  cat("working_directory:", getwd(), "\n")
  cat("command:", paste(commandArgs(), collapse = " "), "\n")
  if (length(inputs)) {
    cat("inputs:\n")
    for (x in inputs) {
      if (file.exists(x)) {
        h <- tryCatch(file_sha256(x), error = function(e) "sha256_error")
        cat(" - ", x, "\t", h, "\n", sep = "")
      } else cat(" - ", x, "\tMISSING\n", sep = "")
    }
  }
  invisible(log_path)
}

finish_log <- function() {
  cat("finished_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
  cat("session_info:\n")
  print(sessionInfo())
  while (sink.number() > 0) sink()
}

safe_mkdir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  safe_mkdir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(as.data.frame(x), path, sep = "\t", quote = FALSE, na = "NA")
  } else {
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
}

file_sha256 <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  if (!requireNamespace("digest", quietly = TRUE)) stop("R package 'digest' is required")
  digest::digest(path, algo = "sha256", file = TRUE)
}

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x <- sub("^gene[-_]", "", x, ignore.case = TRUE)
  x <- sub("^id[-_]", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

as_numeric_safely <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  y[!is.finite(y)] <- NA_real_
  y
}

read_payload <- function(path) {
  stopifnot(file.exists(path))
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    tmp <- tempfile(fileext = ".payload")
    err <- tempfile(fileext = ".stderr")
    status <- system2("gzip", c("-cd", path), stdout = tmp, stderr = err)
    if (file.exists(err)) unlink(err)
    if (!identical(status, 0L) && !identical(status, 0)) stop("gzip failed for ", path)
    on.exit(unlink(tmp), add = TRUE)
    return(readBin(tmp, "raw", n = file.info(tmp)$size))
  }
  readBin(path, "raw", n = file.info(path)$size)
}

payload_type <- function(raw) {
  if (length(raw) >= 4L && identical(raw[1:4], as.raw(c(0x50, 0x4b, 0x03, 0x04)))) return("xlsx")
  if (length(raw) >= 4L && identical(raw[1:4], as.raw(c(0xd0, 0xcf, 0x11, 0xe0)))) return("xls")
  "text"
}

read_geo_table <- function(path) {
  raw <- read_payload(path)
  typ <- payload_type(raw)
  if (typ %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required for ", path)
    ext <- if (typ == "xlsx") ".xlsx" else ".xls"
    tmp <- tempfile(fileext = ext)
    writeBin(raw, tmp)
    on.exit(unlink(tmp), add = TRUE)
    sheets <- readxl::excel_sheets(tmp)
    # GEO tables have one data sheet; use the first sheet and retain all rows.
    tbl <- readxl::read_excel(tmp, sheet = sheets[1], col_names = FALSE, .name_repair = "minimal")
    return(as.data.frame(tbl, check.names = FALSE, stringsAsFactors = FALSE))
  }
  txt <- rawToChar(raw)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  first <- lines[which(nzchar(trimws(lines)))[1]]
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(",", first, fixed = TRUE)) "," else "\t"
  if (requireNamespace("data.table", quietly = TRUE)) {
    tmp <- tempfile(fileext = ".txt")
    writeBin(raw, tmp)
    on.exit(unlink(tmp), add = TRUE)
    return(as.data.frame(data.table::fread(tmp, sep = sep, header = FALSE, fill = TRUE,
                                           quote = "\"", encoding = "UTF-8", data.table = FALSE),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  read.delim(text = txt, sep = sep, header = FALSE, fill = TRUE, quote = "\"",
             check.names = FALSE, stringsAsFactors = FALSE)
}

find_header_row <- function(tbl) {
  # Prefer a true gene/ID header over a title row such as "TPM of all DEGs".
  strong <- apply(tbl, 1, function(r) {
    v <- tolower(trimws(as.character(r)))
    any(grepl("gene[ ._-]*id|gene[ ._-]*name|^id$", v))
  })
  hit <- which(strong)
  if (length(hit)) return(hit[1])
  weak <- apply(tbl, 1, function(r) {
    v <- tolower(trimws(as.character(r)))
    any(grepl("tpm|fpkm|rpkm|count", v))
  })
  hit <- which(weak)
  if (length(hit)) hit[1] else 1L
}

table_with_header <- function(tbl) {
  hrow <- find_header_row(tbl)
  h <- trimws(as.character(unlist(tbl[hrow, , drop = TRUE])))
  h[is.na(h) | !nzchar(h)] <- paste0("V", which(is.na(h) | !nzchar(h)))
  h <- make.unique(h, sep = "__dup")
  body <- if (hrow < nrow(tbl)) tbl[(hrow + 1L):nrow(tbl), , drop = FALSE] else tbl[0, , drop = FALSE]
  colnames(body) <- h
  rownames(body) <- NULL
  attr(body, "header_row") <- hrow
  body
}

map_rows_to_candidates <- function(tbl, candidates, bridge = NULL) {
  candidate_core <- normalize_id(candidates$sequence_id)
  values <- lapply(tbl, normalize_id)
  n <- nrow(tbl)
  mapped <- rep(NA_character_, n)
  method <- rep("not_mapped", n)
  # Direct ID/name matching. This captures AOL_s... and gene-AOL_s... fields.
  for (j in seq_along(values)) {
    m <- match(values[[j]], candidate_core)
    take <- is.na(mapped) & !is.na(m)
    mapped[take] <- candidates$sequence_id[m[take]]
    method[take] <- ifelse(grepl("gene|id|symbol|name", colnames(tbl)[j], ignore.case = TRUE),
                           "direct_id_column", "direct_other_column")
  }
  # Some Majorbio tables retain only internal gene IDs (gene178/MSTRG...).
  if (!is.null(bridge) && nrow(bridge)) {
    id_cols <- grep("gene.?id|^id$", colnames(tbl), ignore.case = TRUE, value = TRUE)
    for (idcol in id_cols) {
      b <- match(normalize_id(tbl[[idcol]]), normalize_id(bridge$raw_gene_id))
      take <- is.na(mapped) & !is.na(b)
      mapped[take] <- bridge$candidate_id[b[take]]
      method[take] <- "gene_id_bridge"
    }
  }
  data.frame(candidate_id = mapped, mapping_method = method, stringsAsFactors = FALSE)
}

make_gene_bridge <- function(parsed_tables, candidates) {
  out <- list()
  for (x in parsed_tables) {
    tbl <- x$table
    m <- map_rows_to_candidates(tbl, candidates, bridge = NULL)
    idcols <- grep("gene.?id|^id$", colnames(tbl), ignore.case = TRUE, value = TRUE)
    namecols <- grep("gene.?name|symbol|^id$", colnames(tbl), ignore.case = TRUE, value = TRUE)
    if (!length(idcols) || !length(namecols)) next
    idv <- as.character(tbl[[idcols[1]]])
    for (nc in namecols) {
      nm <- as.character(tbl[[nc]])
      ok <- !is.na(m$candidate_id) & nzchar(idv) & !is.na(nm) & nzchar(nm)
      if (any(ok)) out[[length(out) + 1L]] <- data.frame(raw_gene_id = idv[ok], candidate_id = m$candidate_id[ok])
    }
  }
  if (!length(out)) return(data.frame(raw_gene_id = character(), candidate_id = character()))
  z <- unique(do.call(rbind, out))
  z[!duplicated(normalize_id(z$raw_gene_id)), , drop = FALSE]
}

parse_time_hours <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  # Prefer an explicit hour token so mutant names such as ham4/chk1/Hal4b
  # cannot be mistaken for the sampling time.
  has_h <- grepl("(?:^|[_-])([0-9]+)h(?:[_-]|$)", x, ignore.case = TRUE, perl = TRUE)
  out[has_h] <- suppressWarnings(as.numeric(sub(".*(?:^|[_-])([0-9]+)h(?:[_-]|$).*", "\\1", x[has_h], ignore.case = TRUE, perl = TRUE)))
  pending <- which(!has_h)
  if (length(pending)) {
    # For no-suffix GEO titles, the time is the first numeric token after the
    # genotype token (e.g. WT_0_1, R7A_12_hypha1, Q405_24_3).
    for (i in pending) {
      parts <- unlist(strsplit(x[i], "[_-]", fixed = FALSE))
      nums <- suppressWarnings(as.numeric(parts))
      j <- which(is.finite(nums))[1]
      if (length(j) && is.finite(j)) {
        # A leading mutant number (Q405, ham4, chk1, Hal4b) is not a time.
        if (j == 1L && length(nums) > 1L) j <- which(is.finite(nums[-1]))[1] + 1L
        if (length(j) && is.finite(j)) out[i] <- nums[j]
      }
    }
  }
  out
}

scale01 <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(rep(0, length(x)))
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || hi == lo) return(ifelse(is.finite(x), 1, 0))
  pmax(0, pmin(1, (x - lo) / (hi - lo)))
}

session_versions <- function() {
  pkgs <- c("readxl", "data.table", "Biostrings", "DESeq2", "ggplot2", "pheatmap", "jsonlite")
  z <- vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, character(1))
  data.frame(package = names(z), version = unname(z), stringsAsFactors = FALSE)
}
