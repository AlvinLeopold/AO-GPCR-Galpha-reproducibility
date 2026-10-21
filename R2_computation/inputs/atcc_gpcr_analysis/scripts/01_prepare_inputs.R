#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(Biostrings)
  library(digest)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

## This script is intentionally self-contained: it reads the user-selected
## workbook, extracts exactly those 93 sequences from ExRe, and freezes all
## input manifests used by downstream steps.
step <- "01_prepare_inputs"
input_paths <- c(
  file.path(REPO, "Aol_GPCR_determined.xlsx"),
  file.path(REPO, "reference", "Arthrobotrys_oligospora.fasta"),
  file.path(REPO, "geo_ao_transcriptomes/catalog/sample_metadata.tsv"),
  file.path(REPO, "geo_ao_transcriptomes/catalog/series_summary.tsv"),
  file.path(PROJ, "config/series_config.tsv")
)
start_log(step, input_paths)

workbook <- file.path(REPO, "Aol_GPCR_determined.xlsx")
exre_fasta <- file.path(REPO, "reference", "Arthrobotrys_oligospora.fasta")
catalog_dir <- file.path(REPO, "geo_ao_transcriptomes", "catalog")
cfg_path <- file.path(PROJ, "config", "series_config.tsv")

candidate_tbl <- as.data.frame(readxl::read_excel(workbook, sheet = "Sheet1", .name_repair = "unique"), stringsAsFactors = FALSE)
required_cols <- c("Sequence_ID", "Species", "TM_Count", "Length", "GPCR_Domain_Count", "GPCR_Domains")
if (!all(required_cols %in% colnames(candidate_tbl))) stop("Workbook columns missing: ", paste(setdiff(required_cols, colnames(candidate_tbl)), collapse = ", "))
candidate_tbl <- candidate_tbl[, required_cols]
colnames(candidate_tbl) <- c("sequence_id", "species", "tm_count", "length_workbook", "gpcr_domain_count", "gpcr_domains")
candidate_tbl$sequence_id <- trimws(as.character(candidate_tbl$sequence_id))
if (nrow(candidate_tbl) != 93L) stop("Expected 93 selected candidates, found ", nrow(candidate_tbl))
if (anyDuplicated(normalize_id(candidate_tbl$sequence_id))) stop("Duplicate candidate IDs after version normalization")
candidate_tbl$candidate_order <- seq_len(nrow(candidate_tbl))
candidate_tbl$candidate_core_id <- normalize_id(candidate_tbl$sequence_id)

seqs <- Biostrings::readAAStringSet(exre_fasta)
seq_headers <- names(seqs)
seq_ids <- sub("\\|.*$", "", seq_headers)
seq_core <- normalize_id(seq_ids)
if (anyDuplicated(seq_core)) stop("Duplicate ExRe sequence IDs after version normalization")
idx <- match(candidate_tbl$candidate_core_id, seq_core)
candidate_tbl$sequence_found_in_exre <- !is.na(idx)
if (any(!candidate_tbl$sequence_found_in_exre)) {
  stop("Candidates missing from ExRe FASTA: ", paste(candidate_tbl$sequence_id[!candidate_tbl$sequence_found_in_exre], collapse = ", "))
}
candidate_tbl$exre_id <- seq_ids[idx]
candidate_tbl$exre_header <- seq_headers[idx]
candidate_tbl$sequence_length <- as.integer(width(seqs[idx]))
candidate_tbl$length_delta <- candidate_tbl$sequence_length - as.numeric(candidate_tbl$length_workbook)
candidate_tbl$sequence_sha256 <- vapply(seqs[idx], function(s) digest(toupper(as.character(s)), algo = "sha256", serialize = FALSE), character(1))
candidate_tbl$sequence_status <- ifelse(abs(candidate_tbl$length_delta) == 0, "length_exact", "length_mismatch_check")

write_tsv(candidate_tbl, file.path(PROJ, "metadata", "candidates.tsv"))

## Write the exact 93-protein FASTA subset with a stable, auditable header.
fasta_lines <- unlist(Map(function(i, s) {
  c(paste0(">", candidate_tbl$sequence_id[i], "|", candidate_tbl$species[i],
           "|TM=", candidate_tbl$tm_count[i], "|GPCR_domains=", candidate_tbl$gpcr_domain_count[i],
           "|sha256=", candidate_tbl$sequence_sha256[i]),
    strsplit(paste(as.character(s), collapse = ""), "(?<=.{80})", perl = TRUE)[[1]])
}, seq_len(nrow(candidate_tbl)), seqs[idx]))
writeLines(fasta_lines, file.path(PROJ, "metadata", "candidate_sequences.fasta"), useBytes = TRUE)

cfg <- fread(cfg_path, sep = "\t", na.strings = c("", "NA"), encoding = "UTF-8")
source_paths <- normalizePath(file.path(PROJ, cfg$file), winslash = "/", mustWork = FALSE)
cfg$file_exists <- file.exists(source_paths)
cfg$file_size_bytes <- ifelse(cfg$file_exists, file.info(source_paths)$size, NA_real_)
cfg$sha256 <- vapply(source_paths, file_sha256, character(1))
cfg$payload_type <- vapply(source_paths, function(p) {
  if (!file.exists(p)) return("missing")
  raw <- read_payload(p)
  payload_type(raw)
}, character(1))
if (any(!cfg$file_exists)) stop("Configured expression files missing: ", paste(cfg$source_id[!cfg$file_exists], collapse = ", "))
write_tsv(cfg, file.path(PROJ, "metadata", "source_file_manifest.tsv"))

sample_meta <- fread(file.path(catalog_dir, "sample_metadata.tsv"), sep = "\t", na.strings = c("", "NA"), encoding = "UTF-8")
series_ids <- unique(cfg$gse)
sample_subset <- sample_meta[gse %in% series_ids]
sample_subset$catalog_source <- "geo_ao_transcriptomes/catalog/sample_metadata.tsv"
write_tsv(sample_subset, file.path(PROJ, "metadata", "sample_metadata_atcc_series.tsv"))

series_summary <- fread(file.path(catalog_dir, "series_summary.tsv"), sep = "\t", na.strings = c("", "NA"), encoding = "UTF-8")
write_tsv(series_summary[gse %in% series_ids], file.path(PROJ, "metadata", "series_summary_atcc.tsv"))

manifest_inputs <- data.frame(
  input_role = c("candidate_workbook", "exre_species_fasta", "geo_sample_metadata", "geo_series_summary", "series_config"),
  path = c("Aol_GPCR_determined.xlsx", "reference/Arthrobotrys_oligospora.fasta",
           "geo_ao_transcriptomes/catalog/sample_metadata.tsv",
           "geo_ao_transcriptomes/catalog/series_summary.tsv",
           "atcc_gpcr_analysis/config/series_config.tsv"),
  exists = file.exists(input_paths),
  size_bytes = ifelse(file.exists(input_paths), file.info(input_paths)$size, NA_real_),
  sha256 = vapply(input_paths, file_sha256, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(manifest_inputs, file.path(PROJ, "checksums", "input_sha256.tsv"))
write_tsv(session_versions(), file.path(PROJ, "metadata", "software_versions.tsv"))

cat("candidate_count:", nrow(candidate_tbl), "\n")
cat("candidate_sequences_found:", sum(candidate_tbl$sequence_found_in_exre), "\n")
cat("length_exact:", sum(candidate_tbl$sequence_status == "length_exact"), "\n")
cat("configured_sources:", nrow(cfg), "\n")
cat("configured_GSEs:", paste(series_ids, collapse = ","), "\n")
cat("sample_count_in_selected_GSEs:", nrow(sample_subset), "\n")
finish_log()
