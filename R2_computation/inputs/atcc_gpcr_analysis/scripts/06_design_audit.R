#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "06_design_audit"
input_paths <- c(
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv"),
  file.path(PROJ, "results", "expression", "candidate_expression_normalized_long.tsv"),
  file.path(PROJ, "results", "expression", "parsed_matrices.rds"),
  file.path(PROJ, "metadata", "source_file_manifest.tsv"),
  file.path(PROJ, "config", "series_config.tsv")
)
start_log(step, input_paths)

dir.create(file.path(PROJ, "results", "design"), recursive = TRUE, showWarnings = FALSE)
candidates <- fread(input_paths[1], sep = "\t", encoding = "UTF-8")
samples <- fread(input_paths[2], sep = "\t", encoding = "UTF-8")
expr <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
parsed <- readRDS(input_paths[4])
source_manifest <- fread(input_paths[5], sep = "\t", encoding = "UTF-8")
series_cfg <- fread(input_paths[6], sep = "\t", encoding = "UTF-8")

if (!all(c("matrix_id", "gse", "gsm", "matrix_column", "genotype", "time_h", "replicate") %in% names(samples))) {
  stop("candidate_sample_metadata.tsv is missing required design columns")
}

safe_paste <- function(x) paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))), collapse = ";")
fmt_time <- function(x) {
  x <- sort(unique(x[is.finite(x)]))
  if (!length(x)) return("Day5_or_unknown")
  paste0(x, "h", collapse = ";")
}
group_key <- function(genotype, time_h) paste(genotype, ifelse(is.finite(time_h), paste0(time_h, "h"), "Day5"), sep = "__")

## Build a design manifest and a fully explicit contrast manifest. Each
## contrast is checked against the sample counts before any model is fitted.
design_rows <- list()
contrast_rows <- list()
mapping_rows <- list()

append_contrast <- function(mid, gse, unit, contrast_type, mutant = NA_character_, time_h = NA_real_,
                            group_a = NA_character_, group_b = NA_character_, note = "") {
  sm <- samples[matrix_id == mid]
  n_a <- if (nzchar(group_a)) sum(group_key(sm$genotype, sm$time_h) == group_a) else 0L
  n_b <- if (nzchar(group_b)) sum(group_key(sm$genotype, sm$time_h) == group_b) else 0L
  has_missing_time <- any(is.na(sm$time_h) & sm$treatment != "aerial_hyphae_differentiation")
  status <- if (has_missing_time) "not_tested" else if (n_a >= 2L && n_b >= 2L) "testable" else "not_tested"
  reason <- if (status == "testable") "" else if (has_missing_time) "unknown_time_or_mapping" else paste0("replicates_below_2:n_a=", n_a, ",n_b=", n_b)
  contrast_rows[[length(contrast_rows) + 1L]] <<- data.frame(
    matrix_id = mid, gse = gse, unit = unit, contrast_id = paste(mid, contrast_type, mutant, time_h, sep = "::"),
    contrast_type = contrast_type, mutant = mutant, time_h = time_h, group_a = group_a, group_b = group_b,
    n_group_a = n_a, n_group_b = n_b, status = status, not_tested_reason = reason,
    mapping_note = note, stringsAsFactors = FALSE
  )
}

for (mid in names(parsed)) {
  rec <- parsed[[mid]]
  sm <- samples[matrix_id == mid]
  if (!nrow(sm)) stop("No sample metadata for ", mid)
  if (ncol(rec$values) != nrow(sm)) stop("Sample count mismatch for ", mid)
  if (!identical(as.character(rec$columns), as.character(sm$matrix_column))) {
    stop("Matrix columns and sample metadata order mismatch for ", mid)
  }
  g <- sort(unique(sm$genotype))
  t <- sort(unique(sm$time_h[is.finite(sm$time_h)]))
  has_wt <- "WT" %in% g
  has_zero <- 0 %in% t
  nonzero_t <- t[t > 0]
  is_count <- identical(tolower(rec$unit), "count")
  mapping <- safe_paste(sm$mapping_confidence)
  model_formula <- if (!length(t) || (length(t) == 1L && !has_zero)) "~ genotype" else "~ genotype + time_factor + genotype:time_factor"
  matrix_status <- if (any(sm$genotype %in% c("unknown", ""))) "review" else "ready"
  if (mapping == "order_based") matrix_status <- "review_order_mapping"
  design_rows[[length(design_rows) + 1L]] <- data.frame(
    matrix_id = mid, gse = rec$gse, unit = rec$unit, expression_type = if (is_count) "raw_counts" else "processed_expression",
    n_genes = nrow(rec$values), n_samples = nrow(sm), n_gsm = uniqueN(sm$gsm),
    genotypes = safe_paste(g), n_genotypes = length(g), time_points = fmt_time(sm$time_h), n_time_points = length(t),
    treatments = safe_paste(sm$treatment), group_replicates_min = min(sm[, .N, by = .(genotype, time_h)]$N),
    model_formula = model_formula, mapping_confidence = mapping, matrix_status = matrix_status,
    source_sha256 = safe_paste(source_manifest[source_id == mid, sha256]), notes = rec$notes,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(nrow(sm))) {
    mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
      matrix_id = mid, gse = rec$gse, matrix_column = sm$matrix_column[j], gsm = sm$gsm[j],
      sample_title = sm$sample_title[j], genotype = sm$genotype[j], time_h = sm$time_h[j],
      time_label = sm$time_label[j], treatment = sm$treatment[j], replicate = sm$replicate[j],
      mapping_confidence = sm$mapping_confidence[j], matrix_column_present = sm$matrix_column[j] %in% rec$columns,
      gsm_unique_within_matrix = sum(sm$gsm == sm$gsm[j]) == 1L, stringsAsFactors = FALSE
    )
  }

  ## WT induction and mutant induction are defined relative to the same
  ## genotype's baseline. Genotype and interaction effects are generated for
  ## every non-baseline time point; baseline genotype effects are retained.
  if (has_wt && has_zero && length(nonzero_t)) {
    for (tt in nonzero_t) {
      append_contrast(mid, rec$gse, rec$unit, "wt_induction", NA_character_, tt,
                      group_key("WT", tt), group_key("WT", 0), "WT(t)-WT(0)")
      for (mut in setdiff(g, "WT")) {
        append_contrast(mid, rec$gse, rec$unit, "mutant_induction", mut, tt,
                        group_key(mut, tt), group_key(mut, 0), "mutant(t)-mutant(0)")
        append_contrast(mid, rec$gse, rec$unit, "genotype_effect", mut, tt,
                        group_key(mut, tt), group_key("WT", tt), "mutant(t)-WT(t)")
        append_contrast(mid, rec$gse, rec$unit, "genotype_time_interaction", mut, tt,
                        group_key(mut, tt), group_key("WT", tt),
                        "[mutant(t)-mutant(0)]-[WT(t)-WT(0)]")
      }
    }
    for (mut in setdiff(g, "WT")) {
      append_contrast(mid, rec$gse, rec$unit, "genotype_effect", mut, 0,
                      group_key(mut, 0), group_key("WT", 0), "baseline mutant(0)-WT(0)")
    }
  } else if (has_wt && length(setdiff(g, "WT"))) {
    for (mut in setdiff(g, "WT")) {
      append_contrast(mid, rec$gse, rec$unit, "genotype_effect", mut, NA_real_,
                      group_key(mut, NA_real_), group_key("WT", NA_real_), "Day5 mutant-WT")
    }
  }
}

design_manifest <- rbindlist(design_rows, fill = TRUE)
contrast_manifest <- rbindlist(contrast_rows, fill = TRUE)
mapping_audit <- rbindlist(mapping_rows, fill = TRUE)
write_tsv(design_manifest, file.path(PROJ, "results", "design", "series_design_manifest.tsv"))
write_tsv(contrast_manifest, file.path(PROJ, "results", "design", "contrast_manifest.tsv"))
write_tsv(mapping_audit, file.path(PROJ, "results", "design", "sample_mapping_audit.tsv"))
write_tsv(contrast_manifest[status == "not_tested"], file.path(PROJ, "results", "design", "contrast_manifest.not_tested.tsv"))

## Explicitly audit possible duplicate GEO submissions at the candidate level.
## The two records have different GSMs but the same design; matching is by
## genotype, time and replicate, never by GSM accession.
dup_a <- expr[gse == "GSE215079" & unit != "count", .(candidate_id, genotype, time_h, replicate, value_a = value)]
dup_b <- expr[gse == "GSE253759" & unit != "count", .(candidate_id, genotype, time_h, replicate, value_b = value)]
dup <- merge(dup_a, dup_b, by = c("candidate_id", "genotype", "time_h", "replicate"), all = TRUE)
dup[, both_finite := is.finite(value_a) & is.finite(value_b)]
dup[, abs_diff := ifelse(both_finite, abs(value_a - value_b), NA_real_)]
dup[, rel_diff := ifelse(both_finite, abs(value_a - value_b) / pmax(abs(value_a), abs(value_b), 1e-12), NA_real_)]
dup_summary <- dup[, {
  ok <- both_finite
  aa <- value_a[ok]; bb <- value_b[ok]
  cor_s <- if (sum(ok) >= 3L && length(unique(aa)) > 1L && length(unique(bb)) > 1L) suppressWarnings(cor(aa, bb, method = "spearman")) else NA_real_
  data.frame(n_pairs = sum(ok), exact_equal_fraction = if (sum(ok)) mean(abs_diff[ok] <= 1e-10) else NA_real_,
             max_abs_diff = if (sum(ok)) max(abs_diff[ok]) else NA_real_, median_abs_diff = if (sum(ok)) median(abs_diff[ok]) else NA_real_,
             median_relative_diff = if (sum(ok)) median(rel_diff[ok]) else NA_real_, spearman = cor_s,
             stringsAsFactors = FALSE)
}, by = candidate_id]
global_dup <- dup[, {
  ok <- both_finite; aa <- value_a[ok]; bb <- value_b[ok]
  data.frame(candidate_id = "ALL_CANDIDATES", n_pairs = sum(ok),
             exact_equal_fraction = if (sum(ok)) mean(abs_diff[ok] <= 1e-10) else NA_real_,
             max_abs_diff = if (sum(ok)) max(abs_diff[ok]) else NA_real_, median_abs_diff = if (sum(ok)) median(abs_diff[ok]) else NA_real_,
             median_relative_diff = if (sum(ok)) median(rel_diff[ok]) else NA_real_,
             spearman = if (sum(ok) >= 3L) suppressWarnings(cor(aa, bb, method = "spearman")) else NA_real_, stringsAsFactors = FALSE)
}]
dup_summary <- rbindlist(list(global_dup, dup_summary), fill = TRUE)
write_tsv(dup_summary, file.path(PROJ, "results", "design", "duplicate_audit.tsv"))

## A GSM overlap audit catches accidental reuse across otherwise distinct
## matrices. It should be zero except for intentional count/processed pairs.
gsm_map <- unique(samples[, .(matrix_id, gse, gsm)])
pair_rows <- list(); mids <- unique(gsm_map$matrix_id)
if (length(mids) >= 2L) for (i in seq_len(length(mids) - 1L)) for (j in (i + 1L):length(mids)) {
  a <- gsm_map[matrix_id == mids[i], gsm]; b <- gsm_map[matrix_id == mids[j], gsm]
  pair_rows[[length(pair_rows) + 1L]] <- data.frame(matrix_a = mids[i], matrix_b = mids[j],
    gse_a = gsm_map[matrix_id == mids[i], gse][1], gse_b = gsm_map[matrix_id == mids[j], gse][1],
    shared_gsm_n = length(intersect(a, b)), shared_gsm = paste(intersect(a, b), collapse = ";"), stringsAsFactors = FALSE)
}
write_tsv(rbindlist(pair_rows, fill = TRUE), file.path(PROJ, "results", "design", "gsm_overlap_audit.tsv"))

## Checksums for the complete Phase 0 artifact set.
out_files <- list.files(file.path(PROJ, "results", "design"), full.names = TRUE)
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE),
          file.path(PROJ, "checksums", "step06_design_sha256.tsv"))
cat("matrices_audited:", nrow(design_manifest), "\n")
cat("contrast_total:", nrow(contrast_manifest), "\n")
cat("contrast_testable:", sum(contrast_manifest$status == "testable"), "\n")
cat("contrast_not_tested:", sum(contrast_manifest$status == "not_tested"), "\n")
cat("duplicate_candidate_rows:", nrow(dup_summary), "\n")
cat("duplicate_global_exact_fraction:", global_dup$exact_equal_fraction, "\n")
cat("gsm_overlap_pairs_with_shared_gsm:", sum(rbindlist(pair_rows, fill = TRUE)$shared_gsm_n > 0), "\n")
finish_log()
