#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "04_integrate_results"
input_paths <- c(
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "results", "expression", "candidate_expression_normalized_long.tsv"),
  file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv"),
  file.path(PROJ, "results", "effects", "candidate_effects.tsv"),
  file.path(PROJ, "results", "effects", "baseline_expression.tsv"),
  file.path(PROJ, "results", "deseq2", "deseq2_candidate_results.tsv"),
  file.path(PROJ, "config", "analysis_parameters.tsv")
)
start_log(step, input_paths)

candidates <- fread(input_paths[1], sep = "\t", encoding = "UTF-8")
expr <- fread(input_paths[2], sep = "\t", encoding = "UTF-8")
samples <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
effects_path <- input_paths[4]
baseline_path <- input_paths[5]
deseq_path <- input_paths[6]
effects <- if (file.exists(effects_path)) fread(effects_path, sep = "\t", encoding = "UTF-8") else data.table()
baseline <- if (file.exists(baseline_path)) fread(baseline_path, sep = "\t", encoding = "UTF-8") else data.table()
deseq <- if (file.exists(deseq_path)) fread(deseq_path, sep = "\t", encoding = "UTF-8") else data.table()
params <- fread(input_paths[7], sep = "\t", encoding = "UTF-8")
param <- function(n, default = NA_character_) {
  x <- params[parameter == n, value]
  if (!length(x)) default else x[1]
}
fdr_cut <- as.numeric(param("effect_fdr_cutoff", "0.05"))
lfc_cut <- as.numeric(param("effect_abs_log2fc_cutoff", "1"))
duplicate_gse <- param("potential_duplicate_gse", "GSE253759")

safe_mean <- function(x) if (!length(x) || all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (!length(x) || all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
consistency <- function(pos, neg) {
  n <- pos + neg
  if (!is.finite(n) || n == 0) return(NA_real_)
  max(pos, neg) / n
}
direction <- function(pos, neg) {
  if ((pos + neg) == 0) return("none")
  if (pos > neg) "up" else if (neg > pos) "down" else "mixed"
}
clean_num <- function(x) { x <- as.numeric(x); x[!is.finite(x)] <- 0; x }

## One row per matrix/GSE describes what can be tested and whether the matrix
## is a potential duplicate-looking record.
series_matrix <- samples[, .(
  gse = first(gse), unit = paste(sort(unique(unit)), collapse = ";"), n_samples = uniqueN(gsm),
  genotypes = paste(sort(unique(genotype)), collapse = ";"), times = paste(sort(unique(time_label)), collapse = ";"),
  treatments = paste(sort(unique(treatment)), collapse = ";"), source_files = paste(unique(source_file), collapse = ";"),
  mapping_confidence = paste(unique(mapping_confidence), collapse = ";")
), by = matrix_id]
series_matrix[, potential_duplicate_record := gse == duplicate_gse]
series_matrix[, sample_count_note := ifelse(n_samples >= 3, "replicates_present", "small_or_unknown")]
write_tsv(series_matrix[order(gse, matrix_id)], file.path(PROJ, "results", "integration", "series_matrix_summary.tsv"))

## Detection breadth: any positive normalized value in a matrix is retained as
## detection; all-zero and absent rows are distinguished by n_observed.
detection <- expr[, .(
  n_observed = sum(is.finite(value)), n_detected = sum(is.finite(value) & value > 0),
  detection_rate = ifelse(sum(is.finite(value)) > 0, mean(is.finite(value) & value > 0), NA_real_),
  mean_value = safe_mean(value), median_value = safe_median(value), mean_log2 = safe_mean(log2(value + 0.1))
), by = .(candidate_id, gse, matrix_id, unit)]
detection[, detected_matrix := n_detected > 0]
detection[, potential_duplicate_record := gse == duplicate_gse]
write_tsv(detection, file.path(PROJ, "results", "integration", "candidate_detection_by_matrix.tsv"))

det_summary <- detection[, .(
  normalized_matrix_n = uniqueN(matrix_id), normalized_gse_n = uniqueN(gse),
  detected_matrix_n = sum(detected_matrix, na.rm = TRUE), detected_gse_n = uniqueN(gse[detected_matrix]),
  detected_gse_n_conservative = uniqueN(gse[detected_matrix & gse != duplicate_gse]),
  median_detected_rate = safe_median(detection_rate), median_mean_log2 = safe_median(mean_log2),
  baseline_like_detected_gse_n = uniqueN(gse[detected_matrix & unit != "count"])
), by = candidate_id]

## WT induction evidence is summarized both by matrix and by distinct GSE. The
## conservative count excludes the potential duplicate-looking GSE253759.
wt <- effects[contrast_type == "wt_induction" & is.finite(log2fc)]
if (nrow(wt)) {
  wt[, direction := ifelse(log2fc > 0, "up", ifelse(log2fc < 0, "down", "flat"))]
  wt[, supported2 := is.finite(p_adj) & p_adj < fdr_cut & abs(log2fc) >= lfc_cut]
  wt_matrix <- wt[, .(
    wt_induction_matrix_n = uniqueN(matrix_id), wt_induction_supported_matrix_n = uniqueN(matrix_id[supported2]),
    wt_induction_gse_n = uniqueN(gse), wt_induction_supported_gse_n = uniqueN(gse[supported2]),
    wt_induction_supported_gse_n_conservative = uniqueN(gse[supported2 & gse != duplicate_gse]),
    wt_up_n = sum(log2fc > 0), wt_down_n = sum(log2fc < 0), wt_supported_up_n = sum(supported2 & log2fc > 0),
    wt_supported_down_n = sum(supported2 & log2fc < 0), wt_max_abs_log2fc = max(abs(log2fc), na.rm = TRUE),
    wt_median_log2fc = safe_median(log2fc)
  ), by = candidate_id]
  wt_gse <- wt[, .(
    gse_log2fc = safe_mean(log2fc), gse_supported = any(supported2), gse_direction = direction(sum(log2fc > 0), sum(log2fc < 0)),
    gse_potential_duplicate = first(gse == duplicate_gse)
  ), by = .(candidate_id, gse)]
  wt_gse_summary <- wt_gse[, .(
  wt_gse_replicated_n = .N, wt_gse_supported_n = sum(gse_supported), wt_gse_supported_n_conservative = sum(gse_supported & !gse_potential_duplicate),
    wt_gse_up_n = sum(gse_direction == "up"), wt_gse_down_n = sum(gse_direction == "down"),
    wt_gse_direction_consistency = consistency(sum(gse_direction == "up"), sum(gse_direction == "down")),
    wt_gse_direction = direction(sum(gse_direction == "up"), sum(gse_direction == "down")),
    wt_gse_up_n_conservative = sum(gse_direction == "up" & !gse_potential_duplicate),
    wt_gse_down_n_conservative = sum(gse_direction == "down" & !gse_potential_duplicate),
    wt_gse_direction_consistency_conservative = consistency(sum(gse_direction == "up" & !gse_potential_duplicate), sum(gse_direction == "down" & !gse_potential_duplicate)),
    wt_gse_direction_conservative = direction(sum(gse_direction == "up" & !gse_potential_duplicate), sum(gse_direction == "down" & !gse_potential_duplicate))
  ), by = candidate_id]
} else {
  wt_matrix <- data.table(); wt_gse_summary <- data.table()
}

gen <- effects[contrast_type == "genotype_effect" & is.finite(log2fc)]
if (nrow(gen)) {
  gen[, supported2 := is.finite(p_adj) & p_adj < fdr_cut & abs(log2fc) >= lfc_cut]
  gen_summary <- gen[, .(
    genotype_effect_matrix_n = uniqueN(matrix_id), genotype_effect_supported_matrix_n = uniqueN(matrix_id[supported2]),
    genotype_effect_gse_n = uniqueN(gse), genotype_effect_supported_gse_n = uniqueN(gse[supported2]),
    genotype_effect_supported_gse_n_conservative = uniqueN(gse[supported2 & gse != duplicate_gse]),
    genotype_effect_up_n = sum(supported2 & log2fc > 0), genotype_effect_down_n = sum(supported2 & log2fc < 0),
    genotype_effect_max_abs_log2fc = max(abs(log2fc), na.rm = TRUE)
  ), by = candidate_id]
} else gen_summary <- data.table()

int <- effects[contrast_type == "genotype_time_interaction" & is.finite(log2fc)]
if (nrow(int)) {
  int[, supported2 := is.finite(p_adj) & p_adj < fdr_cut & abs(log2fc) >= lfc_cut]
  int_summary <- int[, .(
    interaction_matrix_n = uniqueN(matrix_id), interaction_supported_matrix_n = uniqueN(matrix_id[supported2]),
    interaction_gse_n = uniqueN(gse), interaction_supported_gse_n = uniqueN(gse[supported2]),
    interaction_supported_gse_n_conservative = uniqueN(gse[supported2 & gse != duplicate_gse]),
    interaction_max_abs_log2fc = max(abs(log2fc), na.rm = TRUE)
  ), by = candidate_id]
} else int_summary <- data.table()

if (nrow(deseq)) {
  deseq[, supported2 := is.finite(padj) & padj < fdr_cut & abs(log2FoldChange) >= lfc_cut]
  de_summary <- deseq[, .(
    deseq_result_n = .N, deseq_supported_n = sum(supported2), deseq_gse_n = uniqueN(gse),
    deseq_supported_gse_n = uniqueN(gse[supported2]), deseq_supported_gse_n_conservative = uniqueN(gse[supported2 & gse != duplicate_gse]),
    deseq_interaction_supported_gse_n = uniqueN(gse[supported2 & contrast_type == "DESeq2_interaction"]),
    deseq_max_abs_log2fc = max(abs(log2FoldChange), na.rm = TRUE)
  ), by = candidate_id]
} else de_summary <- data.table()

## Evidence by GSE: a compact matrix for manuscript/table review. Missing data
## are explicit rather than converted to zeros.
gse_order <- sort(unique(c(samples$gse)))
gse_rows <- list()
for (cid in candidates$sequence_id) {
  for (g in gse_order) {
    dd <- detection[candidate_id == cid & gse == g]
    ww <- wt[candidate_id == cid & gse == g]
    gg <- gen[candidate_id == cid & gse == g]
    ii <- int[candidate_id == cid & gse == g]
    de <- deseq[candidate_id == cid & gse == g]
    gse_rows[[length(gse_rows) + 1L]] <- data.frame(candidate_id = cid, gse = g,
      baseline_detected = any(dd$detected_matrix[dd$unit != "count"], na.rm = TRUE),
      baseline_matrix_n = uniqueN(dd$matrix_id[dd$unit != "count"]),
      wt_induction_records = nrow(ww), wt_induction_supported = sum(ww$supported2, na.rm = TRUE),
      wt_induction_direction = if (nrow(ww)) direction(sum(ww$log2fc > 0), sum(ww$log2fc < 0)) else "none",
      wt_induction_median_log2fc = if (nrow(ww)) safe_median(ww$log2fc) else NA_real_,
      genotype_effect_records = nrow(gg), genotype_effect_supported = sum(gg$supported2, na.rm = TRUE),
      interaction_records = nrow(ii), interaction_supported = sum(ii$supported2, na.rm = TRUE),
      deseq_records = nrow(de), deseq_supported = sum(de$supported2, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
}
evidence_by_gse <- rbindlist(gse_rows, fill = TRUE)
write_tsv(evidence_by_gse, file.path(PROJ, "results", "integration", "candidate_evidence_by_gse.tsv"))

summary <- merge(candidates[, .(candidate_id = sequence_id, sequence_id, candidate_order, tm_count, length_workbook,
                               gpcr_domain_count, gpcr_domains, sequence_length, sequence_sha256)], det_summary, by = "candidate_id", all.x = TRUE)
for (x in list(wt_matrix, wt_gse_summary, gen_summary, int_summary, de_summary)) if (nrow(x)) summary <- merge(summary, x, by = "candidate_id", all.x = TRUE)
for (j in names(summary)) if (is.numeric(summary[[j]])) summary[[j]][!is.finite(summary[[j]])] <- 0
summary[, baseline_breadth_score := pmin(1, detected_gse_n_conservative / 5)]
summary[, abundance_score := scale01(median_mean_log2)]
summary[, induction_support_score := pmin(1, wt_induction_supported_gse_n_conservative / 3)]
summary[, induction_consistency_score := fifelse(wt_gse_direction_consistency_conservative > 0, wt_gse_direction_consistency_conservative, 0)]
summary[, genotype_support_score := pmin(1, genotype_effect_supported_gse_n_conservative / 2)]
summary[, interaction_support_score := pmin(1, interaction_supported_gse_n_conservative / 2)]
summary[, deseq_support_score := pmin(1, deseq_supported_gse_n_conservative / 2)]
summary[, priority_score := 0.15 * abundance_score + 0.15 * baseline_breadth_score +
         0.25 * induction_support_score + 0.15 * induction_consistency_score +
         0.15 * genotype_support_score + 0.10 * interaction_support_score + 0.05 * deseq_support_score]
summary[, priority_tier := fifelse(
  (detected_gse_n_conservative >= as.numeric(param("tier1_min_baseline_gse", "3")) &
     wt_induction_supported_gse_n_conservative >= 3 &
     wt_gse_direction_consistency_conservative >= as.numeric(param("tier1_min_induction_consistency", "0.67"))), "Tier 1: multi-GSE supported",
  fifelse(detected_gse_n_conservative >= as.numeric(param("tier2_min_baseline_gse", "2")) &
            ((wt_induction_supported_gse_n_conservative >= 2 & wt_gse_direction_consistency_conservative >= as.numeric(param("tier2_min_induction_consistency", "0.6"))) |
             (deseq_supported_gse_n_conservative >= as.numeric(param("tier1_min_deseq_gse", "2")) & genotype_effect_supported_gse_n_conservative >= 1)),
          "Tier 2: focused follow-up", "Tier 3: exploratory/low support"))]
summary[, recommended_validation := fifelse(priority_tier == "Tier 1: multi-GSE supported",
  "qRT-PCR/targeted expression plus trap-induction phenotype; test receptor localization or ligand response",
  fifelse(priority_tier == "Tier 2: focused follow-up", "confirm expression in an independent ATCC culture and test condition-specific phenotype",
          "retain as exploratory candidate; first verify transcript detection and gene model"))]
setorder(summary, -priority_score, candidate_order)
summary[, priority_rank := seq_len(.N)]
write_tsv(summary, file.path(PROJ, "results", "integration", "candidate_cross_study_summary.tsv"))

## Direction heatmap: one cell per candidate/GSE, WT induction median log2FC.
if (nrow(wt)) {
  hm <- wt[, .(log2fc = safe_median(log2fc)), by = .(candidate_id, gse)]
  hm <- dcast(hm, candidate_id ~ gse, value.var = "log2fc")
  mat <- as.matrix(hm[, -1, with = FALSE]); rownames(mat) <- hm$candidate_id
  mat[!is.finite(mat)] <- NA_real_
  png(file.path(PROJ, "figures", "wt_induction_direction_heatmap.png"), width = 2100, height = 2800, res = 180)
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, na_col = "grey90",
           main = "WT induction response across GEO records (median log2FC)")
  dev.off()
}

top <- summary[1:min(25L, .N)]
top_long <- melt(top, id.vars = "candidate_id", measure.vars = c("abundance_score", "baseline_breadth_score", "induction_support_score",
                                                                    "induction_consistency_score", "genotype_support_score", "interaction_support_score", "deseq_support_score"),
                 variable.name = "evidence_component", value.name = "component_score")
png(file.path(PROJ, "figures", "top25_evidence_components.png"), width = 2200, height = 1500, res = 180)
print(ggplot(top_long, aes(x = reorder(candidate_id, -component_score), y = component_score, fill = evidence_component)) +
        geom_col(position = "dodge") + coord_flip() + theme_bw(base_size = 8) +
        labs(x = NULL, y = "component score", title = "Evidence components for top 25 candidates") + theme(legend.position = "bottom"))
dev.off()

write_tsv(data.frame(path = list.files(file.path(PROJ, "figures"), full.names = TRUE),
                     sha256 = vapply(list.files(file.path(PROJ, "figures"), full.names = TRUE), file_sha256, character(1))),
          file.path(PROJ, "checksums", "step04_figures_sha256.tsv"))
write_tsv(summary[, .(candidate_id, priority_rank, priority_tier, priority_score, detected_gse_n_conservative,
                       wt_induction_supported_gse_n_conservative, genotype_effect_supported_gse_n_conservative,
                       interaction_supported_gse_n_conservative, deseq_supported_gse_n_conservative,
                       wt_gse_direction, wt_gse_direction_consistency, wt_gse_direction_conservative, wt_gse_direction_consistency_conservative, recommended_validation)],
          file.path(PROJ, "results", "integration", "candidate_priority_shortlist.tsv"))

cat("candidate_count:", nrow(summary), "\n")
cat("GSE_count:", length(gse_order), "\n")
cat("tier_counts:\n"); print(summary[, .N, by = priority_tier])
cat("top_candidates:\n"); print(summary[1:min(20L, .N), .(priority_rank, candidate_id, priority_tier, priority_score, detected_gse_n_conservative,
  wt_induction_supported_gse_n_conservative, genotype_effect_supported_gse_n_conservative, deseq_supported_gse_n_conservative, wt_gse_direction_conservative, wt_gse_direction_consistency_conservative)])
finish_log()
