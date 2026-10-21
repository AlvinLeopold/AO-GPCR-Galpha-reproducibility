#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "11_condition_panels"
input_paths <- c(
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "results", "perturbation", "candidate_regulator_summary.tsv"),
  file.path(PROJ, "results", "time_response", "wt_time_consensus.tsv"),
  file.path(PROJ, "results", "network", "candidate_module_consensus.tsv"),
  file.path(PROJ, "results", "models", "candidate_series_model_results.tsv"),
  file.path(PROJ, "config", "deep_analysis_parameters.tsv")
)
start_log(step, input_paths)
candidates <- fread(input_paths[1], sep = "\t", encoding = "UTF-8")
reg <- fread(input_paths[2], sep = "\t", encoding = "UTF-8")
time <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
net <- if (file.exists(input_paths[4])) fread(input_paths[4], sep = "\t", encoding = "UTF-8") else data.table()
models <- fread(input_paths[5], sep = "\t", encoding = "UTF-8")
params <- fread(input_paths[6], sep = "\t", encoding = "UTF-8")
param <- function(n, default) { x <- params[parameter == n, value]; if (length(x)) x[1] else default }
panel_size <- as.integer(param("panel_size", "10"))

safe_mkdir(file.path(PROJ, "results", "panels")); safe_mkdir(file.path(PROJ, "figures", "panels"))
base <- candidates[, .(candidate_id = sequence_id, candidate_order, tm_count, sequence_length, gpcr_domain_count, gpcr_domains)]

## Convert the regulator summary into candidate-level evidence features.
feat <- reg[, .(
  wt_supported = sum(wt_induction_supported_n > 0), interaction_supported = sum(interaction_supported_n > 0),
  preserved = sum(dependency_class == "mutant_induction_preserved"), lost_or_changed = sum(dependency_class == "mutant_induction_lost_or_changed"),
  baseline_shift = sum(dependency_class == "baseline_genotype_shift"), max_wt = if (any(is.finite(wt_induction_max_abs_log2fc))) max(wt_induction_max_abs_log2fc, na.rm = TRUE) else 0,
  max_interaction = if (any(is.finite(interaction_max_abs_log2fc))) max(interaction_max_abs_log2fc, na.rm = TRUE) else 0, regulators = paste(regulator[interaction_supported_n > 0], collapse = ";"),
  pathways = paste(unique(pathway[interaction_supported_n > 0]), collapse = ";")
), by = candidate_id]
for (j in names(feat)) if (is.numeric(feat[[j]])) feat[[j]][!is.finite(feat[[j]])] <- 0

## Raw-count evidence is limited to strong candidate rows and summarized by
## distinct GEO accession to avoid double-counting processed/count pairs.
de_feat <- models[unit == "count" & is.finite(p_adj), .(
  count_gse_supported = uniqueN(gse[p_adj < 0.05 & abs(log2fc) >= 1]), count_max_abs_log2fc = if (any(p_adj < 0.05)) max(abs(log2fc[p_adj < 0.05]), na.rm = TRUE) else 0
), by = candidate_id]
time_feat <- time[, .(time_class = paste(unique(na.omit(consensus_class)), collapse = ";"), time_peak_h = if (any(is.finite(consensus_peak_time_h))) as.numeric(median(consensus_peak_time_h, na.rm = TRUE)) else NA_real_, time_max_abs = if (any(is.finite(consensus_max_abs_log2fc))) as.numeric(max(consensus_max_abs_log2fc, na.rm = TRUE)) else 0), by = candidate_id]
for (z in list(de_feat, time_feat)) if (nrow(z)) for (j in names(z)) if (is.numeric(z[[j]])) z[[j]][!is.finite(z[[j]])] <- 0
all_feat <- Reduce(function(a, b) merge(a, b, by = "candidate_id", all = TRUE), list(base, feat, de_feat, time_feat))
for (j in names(all_feat)) if (is.numeric(all_feat[[j]])) all_feat[[j]][!is.finite(all_feat[[j]])] <- 0
all_feat[is.na(all_feat)] <- 0

## Panel scores are intentionally condition-specific. They prioritize effect
## size, replication, and the relevant perturbation axis; they are not claims
## of receptor function.
panel_defs <- data.table(
  panel_id = c("P1_predation_core", "P2_cAMP_PKA", "P3_MAPK_osmotic", "P4_vacuolar_trafficking", "P5_hyphal_fusion", "P6_nutrition_autophagy_transcription"),
  focus = c("WT nematode induction and time response", "AoAcy/AoPkaC1 dependence", "AoOpy2/AoHog1 dependence", "AoRab7/AoVps35 dependence", "AoHam4/AoHam8/AoHal4b dependence", "AoAtg11/AobZip129/AoAge2 dependence"),
  stringsAsFactors = FALSE
)
score_panel <- function(d, pid) {
  common <- 0.15 * scale01(d$wt_supported) + 0.15 * scale01(d$max_wt) + 0.15 * scale01(d$count_gse_supported) + 0.10 * scale01(d$time_max_abs)
  if (pid == "P1_predation_core") return(common + 0.20 * scale01(d$wt_supported) + 0.15 * scale01(d$time_max_abs) + 0.10 * scale01(d$max_wt))
  if (pid == "P2_cAMP_PKA") return(common + 0.25 * scale01(d$interaction_supported) + 0.15 * scale01(d$preserved) + 0.10 * scale01(d$lost_or_changed))
  if (pid == "P3_MAPK_osmotic") return(common + 0.25 * scale01(d$interaction_supported) + 0.15 * scale01(d$max_interaction) + 0.10 * scale01(d$baseline_shift))
  if (pid == "P4_vacuolar_trafficking") return(common + 0.30 * scale01(d$interaction_supported) + 0.15 * scale01(d$preserved) + 0.10 * scale01(d$max_interaction))
  if (pid == "P5_hyphal_fusion") return(common + 0.25 * scale01(d$interaction_supported) + 0.15 * scale01(d$lost_or_changed) + 0.10 * scale01(d$baseline_shift))
  common + 0.20 * scale01(d$interaction_supported) + 0.15 * scale01(d$baseline_shift) + 0.10 * scale01(d$preserved)
}

panel_rows <- list(); panel_top <- list()
for (i in seq_len(nrow(panel_defs))) {
  pid <- panel_defs$panel_id[i]; score <- score_panel(all_feat, pid)
  d <- copy(all_feat); d[, `:=`(panel_id = pid, panel_focus = panel_defs$focus[i], panel_score = score)]
  d[, panel_rank := frank(-panel_score, ties.method = "first")]
  d[, panel_selected := panel_rank <= panel_size]
  d[, selection_basis := paste(panel_defs$focus[i], "; score combines effect size, replication, and axis-specific perturbation evidence")]
  panel_rows[[length(panel_rows) + 1L]] <- d
  panel_top[[length(panel_top) + 1L]] <- d[panel_selected == TRUE, .(panel_id, panel_focus, panel_rank, candidate_id, panel_score, tm_count, sequence_length, wt_supported, interaction_supported, preserved, lost_or_changed, baseline_shift, count_gse_supported, time_class, time_peak_h, regulators, pathways, selection_basis)]
}
panel_long <- rbindlist(panel_rows, fill = TRUE); panel_selected <- rbindlist(panel_top, fill = TRUE)
setorder(panel_selected, panel_id, panel_rank)
write_tsv(panel_long, file.path(PROJ, "results", "panels", "candidate_panel_scores_long.tsv"))
write_tsv(panel_selected, file.path(PROJ, "results", "panels", "condition_specific_candidate_panels.tsv"))
write_tsv(panel_defs, file.path(PROJ, "results", "panels", "panel_definitions.tsv"))

## A compact overview matrix records panel membership and ranks.
overview <- dcast(panel_selected, candidate_id ~ panel_id, value.var = "panel_rank")
write_tsv(overview, file.path(PROJ, "results", "panels", "candidate_panel_membership.tsv"))

out_files <- list.files(file.path(PROJ, "results", "panels"), full.names = TRUE)
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE), file.path(PROJ, "checksums", "step11_panels_sha256.tsv"))
cat("panel_count:", nrow(panel_defs), "\n"); cat("panel_selected_rows:", nrow(panel_selected), "\n"); cat("unique_candidates_selected:", uniqueN(panel_selected$candidate_id), "\n")
cat("panel_sizes:\n"); print(panel_selected[, .N, by = panel_id])
finish_log()
