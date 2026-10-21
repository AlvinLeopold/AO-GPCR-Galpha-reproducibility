#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "09_perturbation_consensus"
input_paths <- c(
  file.path(PROJ, "results", "models", "candidate_series_model_results.tsv"),
  file.path(PROJ, "results", "design", "duplicate_audit.tsv"),
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "config", "deep_analysis_parameters.tsv")
)
start_log(step, input_paths)

models <- fread(input_paths[1], sep = "\t", encoding = "UTF-8")
candidates <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
params <- fread(input_paths[4], sep = "\t", encoding = "UTF-8")
param <- function(n, default) { x <- params[parameter == n, value]; if (length(x)) x[1] else default }
fdr_cut <- as.numeric(param("fdr_cutoff", "0.05"))
lfc_cut <- as.numeric(param("abs_log2fc_cutoff", "1"))
duplicate_gse <- param("duplicate_gse", "GSE253759")

safe_mkdir(file.path(PROJ, "results", "perturbation"))
safe_mkdir(file.path(PROJ, "figures", "perturbation"))

regulator_map <- data.table(
  matrix_id = c("GSE210886_Acy", "GSE210886_PKA", "GSE212735", "GSE213542", "GSE215079", "GSE245014_Rab7", "GSE245014_Vps35",
                "GSE301423_counts", "GSE301452_counts", "GSE303416_counts", "GSE306156_counts", "GSE322548_counts", "GSE322740_FPKM"),
  regulator = c("AoAcy", "AoPkaC1", "AoAtg11", "AobZip129", "AoAge2", "AoRab7", "AoVps35", "AoHam4_or_AoHam8", "AoChk1", "AoHal4b", "AoOpy2", "AoHog1", "AoMdv1_or_AoVps1"),
  pathway = c("cAMP-adenylyl cyclase", "cAMP-PKA", "nutrition-autophagy", "transcriptional regulation", "Arf-GAP/Age2", "vacuolar trafficking", "vacuolar trafficking", "hyphal fusion", "checkpoint", "hyphal fusion", "osmosensing-MAPK", "osmosensing-MAPK", "mitochondrial fission"),
  stringsAsFactors = FALSE
)
## Split multi-mutant matrices into the corresponding regulator label below.
regulator_for <- function(mid, mutant) {
  if (mid == "GSE301423_counts") return(ifelse(grepl("ham4", mutant, ignore.case = TRUE), "AoHam4", "AoHam8"))
  if (mid == "GSE322740_FPKM") return(ifelse(grepl("Mdv1", mutant, ignore.case = TRUE), "AoMdv1", "AoVps1"))
  regulator_map[matrix_id == mid, regulator][1]
}
pathway_for <- function(reg) {
  if (grepl("Acy|Pka", reg, ignore.case = TRUE)) return("cAMP-PKA")
  if (grepl("Atg11", reg, ignore.case = TRUE)) return("nutrition-autophagy")
  if (grepl("bZip", reg, ignore.case = TRUE)) return("transcriptional regulation")
  if (grepl("Age2", reg, ignore.case = TRUE)) return("Arf-GAP/Age2")
  if (grepl("Rab7|Vps35", reg, ignore.case = TRUE)) return("vacuolar trafficking")
  if (grepl("Ham4|Ham8|Hal4b", reg, ignore.case = TRUE)) return("hyphal fusion")
  if (grepl("Chk1", reg, ignore.case = TRUE)) return("checkpoint")
  if (grepl("Opy2|Hog1", reg, ignore.case = TRUE)) return("osmosensing-MAPK")
  if (grepl("Mdv1|Vps1", reg, ignore.case = TRUE)) return("mitochondrial fission")
  "other"
}

## Prefer raw-count confirmation when both a processed and count matrix are
## available for the same GEO design. The duplicate Age2 GEO record is retained
## in the evidence table but excluded from independent-support counts.
models[, source_priority := ifelse(unit == "count", 2L, 1L)]
models[, regulator := vapply(seq_len(.N), function(i) regulator_for(matrix_id[i], mutant[i]), character(1))]
models[, pathway := vapply(regulator, pathway_for, character(1))]
models[, independent_gse := gse != duplicate_gse]
setorder(models, candidate_id, regulator, gse, contrast_type, mutant, time_h, -source_priority)
best <- models[, .SD[1], by = .(candidate_id, regulator, gse, contrast_type, mutant, time_h)]
best[, supported_strong := is.finite(p_adj) & p_adj < fdr_cut & is.finite(log2fc) & abs(log2fc) >= lfc_cut]
best[, supported_fdr := is.finite(p_adj) & p_adj < fdr_cut]
best[, direction := fifelse(log2fc > 0, "up", fifelse(log2fc < 0, "down", "flat"))]
write_tsv(best[, .(candidate_id, regulator, pathway, matrix_id, gse, unit, method, contrast_type, mutant, time_h, log2fc, se, ci_low, ci_high, p_value, p_adj,
                   supported_fdr, supported_strong, direction, independent_gse, provisional_mapping, mapping_confidence, model_role)],
          file.path(PROJ, "results", "perturbation", "perturbation_evidence_long.tsv"))

max_abs <- function(x) if (any(is.finite(x))) max(abs(x), na.rm = TRUE) else NA_real_
median_num <- function(x) if (any(is.finite(x))) median(x, na.rm = TRUE) else NA_real_
support_direction <- function(d) {
  d <- d[is.finite(log2fc)]
  if (!nrow(d)) return("none")
  pos <- sum(d$supported_strong & d$log2fc > 0, na.rm = TRUE); neg <- sum(d$supported_strong & d$log2fc < 0, na.rm = TRUE)
  if (!pos && !neg) "none" else if (pos > neg) "up" else if (neg > pos) "down" else "mixed"
}

summary_rows <- list()
for (cid in candidates$sequence_id) for (reg in unique(best$regulator)) {
  d <- best[candidate_id == cid & regulator == reg]
  wt <- d[contrast_type == "wt_induction" & independent_gse]
  mu <- d[contrast_type == "mutant_induction" & independent_gse]
  ge <- d[contrast_type == "genotype_effect" & independent_gse]
  it <- d[contrast_type == "genotype_time_interaction" & independent_gse]
  wt_strong <- any(wt$supported_strong, na.rm = TRUE); mu_strong <- any(mu$supported_strong, na.rm = TRUE); it_strong <- any(it$supported_strong, na.rm = TRUE)
  wt_dir <- support_direction(wt); mu_dir <- support_direction(mu); it_dir <- support_direction(it)
  dep <- if (it_strong) "interaction_supported" else if (wt_strong && mu_strong && wt_dir == mu_dir) "mutant_induction_preserved" else if (wt_strong && (!mu_strong || mu_dir != wt_dir)) "mutant_induction_lost_or_changed" else if (any(ge$supported_strong, na.rm = TRUE)) "baseline_genotype_shift" else "weak_or_unresolved"
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(candidate_id = cid, regulator = reg, pathway = pathway_for(reg),
    wt_induction_n = nrow(wt), wt_induction_supported_n = sum(wt$supported_strong, na.rm = TRUE), wt_induction_direction = wt_dir,
    wt_induction_max_abs_log2fc = max_abs(wt$log2fc), wt_induction_median_log2fc = median_num(wt$log2fc),
    mutant_induction_n = nrow(mu), mutant_induction_supported_n = sum(mu$supported_strong, na.rm = TRUE), mutant_induction_direction = mu_dir,
    mutant_induction_max_abs_log2fc = max_abs(mu$log2fc), mutant_induction_median_log2fc = median_num(mu$log2fc),
    genotype_effect_n = nrow(ge), genotype_effect_supported_n = sum(ge$supported_strong, na.rm = TRUE), genotype_effect_max_abs_log2fc = max_abs(ge$log2fc),
    interaction_n = nrow(it), interaction_supported_n = sum(it$supported_strong, na.rm = TRUE), interaction_direction = it_dir,
    interaction_max_abs_log2fc = max_abs(it$log2fc), interaction_median_log2fc = median_num(it$log2fc),
    independent_gse_n = uniqueN(d$gse[d$independent_gse]), dependency_class = dep,
    provisional_mapping = any(d$provisional_mapping, na.rm = TRUE), stringsAsFactors = FALSE)
}
reg_summary <- rbindlist(summary_rows, fill = TRUE)
write_tsv(reg_summary, file.path(PROJ, "results", "perturbation", "candidate_regulator_summary.tsv"))

## Pairwise perturbation comparisons use the same candidate and time point. A
## pair is called concordant only when the two effect directions agree; absence
## of a result remains explicit and is never treated as zero.
pair_defs <- list(
  c("AoAcy", "AoPkaC1", "cAMP-PKA"), c("AoRab7", "AoVps35", "vacuolar-trafficking"), c("AoHam4", "AoHam8", "hyphal-fusion")
)
pair_rows <- list()
for (pd in pair_defs) {
  a <- best[regulator == pd[1] & contrast_type %in% c("wt_induction", "genotype_time_interaction") & independent_gse,
            .(candidate_id, contrast_type, time_h, log2fc_a = log2fc, supported_a = supported_strong, p_adj_a = p_adj)]
  b <- best[regulator == pd[2] & contrast_type %in% c("wt_induction", "genotype_time_interaction") & independent_gse,
            .(candidate_id, contrast_type, time_h, log2fc_b = log2fc, supported_b = supported_strong, p_adj_b = p_adj)]
  z <- merge(a, b, by = c("candidate_id", "contrast_type", "time_h"), all = TRUE)
  z[, `:=`(regulator_a = pd[1], regulator_b = pd[2], pathway = pd[3], direction_agree = ifelse(is.finite(log2fc_a) & is.finite(log2fc_b), sign(log2fc_a) == sign(log2fc_b), NA),
           joint_supported = supported_a & supported_b, specificity = fifelse(supported_a & !supported_b, paste0(pd[1], "_specific"), fifelse(supported_b & !supported_a, paste0(pd[2], "_specific"), fifelse(supported_a & supported_b, "shared", "unresolved"))))]
  pair_rows[[length(pair_rows) + 1L]] <- z
}
pairwise <- rbindlist(pair_rows, fill = TRUE)
write_tsv(pairwise, file.path(PROJ, "results", "perturbation", "paired_perturbation_comparison.tsv"))

## Candidate-level pathway fingerprint, useful for selecting condition-specific
## validation panels without collapsing all biological questions into one score.
fingerprint <- reg_summary[, .(
  regulator_n = uniqueN(regulator), interaction_supported_regulator_n = sum(interaction_supported_n > 0),
  wt_supported_regulator_n = sum(wt_induction_supported_n > 0), mutant_preserved_n = sum(dependency_class == "mutant_induction_preserved"),
  mutant_lost_or_changed_n = sum(dependency_class == "mutant_induction_lost_or_changed"), baseline_shift_n = sum(dependency_class == "baseline_genotype_shift"),
  interaction_regulators = paste(regulator[interaction_supported_n > 0], collapse = ";"), preserved_regulators = paste(regulator[dependency_class == "mutant_induction_preserved"], collapse = ";"),
  lost_or_changed_regulators = paste(regulator[dependency_class == "mutant_induction_lost_or_changed"], collapse = ";")
), by = candidate_id]
write_tsv(merge(candidates[, .(candidate_id = sequence_id, candidate_order)], fingerprint, by = "candidate_id", all.x = TRUE),
          file.path(PROJ, "results", "perturbation", "candidate_pathway_fingerprint.tsv"))

## Heatmap of interaction magnitudes for the 93-candidate set.
hm <- reg_summary[, .(candidate_id, regulator, value = interaction_max_abs_log2fc)]
hm <- dcast(hm, candidate_id ~ regulator, value.var = "value")
mat <- as.matrix(hm[, -1, with = FALSE]); rownames(mat) <- hm$candidate_id
png(file.path(PROJ, "figures", "perturbation", "interaction_magnitude_heatmap.png"), width = 2200, height = 2600, res = 180)
if (nrow(mat)) {
  ## Missing interaction estimates are shown as zero only in this display;
  ## the tabular outputs preserve NA and the distinction from a tested zero.
  mat_plot <- mat; mat_plot[!is.finite(mat_plot)] <- 0
  suppressWarnings(pheatmap::pheatmap(mat_plot, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, na_col = "grey90", main = "GPCR interaction magnitude by perturbation"))
}
dev.off()

out_files <- c(list.files(file.path(PROJ, "results", "perturbation"), full.names = TRUE), list.files(file.path(PROJ, "figures", "perturbation"), full.names = TRUE))
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE),
          file.path(PROJ, "checksums", "step09_perturbation_sha256.tsv"))
cat("best_evidence_rows:", nrow(best), "\n")
cat("candidate_regulator_rows:", nrow(reg_summary), "\n")
cat("interaction_supported_rows:", sum(reg_summary$interaction_supported_n > 0), "\n")
cat("dependency_class_counts:\n"); print(reg_summary[, .N, by = dependency_class])
cat("pairwise_rows:", nrow(pairwise), "\n")
finish_log()
