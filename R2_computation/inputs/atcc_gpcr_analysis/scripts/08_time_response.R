#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "08_time_response"
input_paths <- c(
  file.path(PROJ, "results", "expression", "candidate_expression_normalized_long.tsv"),
  file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv"),
  file.path(PROJ, "results", "models", "candidate_series_model_results.tsv"),
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "config", "deep_analysis_parameters.tsv")
)
start_log(step, input_paths)

expr <- fread(input_paths[1], sep = "\t", encoding = "UTF-8")
samples <- fread(input_paths[2], sep = "\t", encoding = "UTF-8")
models <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
candidates <- fread(input_paths[4], sep = "\t", encoding = "UTF-8")
params <- fread(input_paths[5], sep = "\t", encoding = "UTF-8")
param <- function(n, default) { x <- params[parameter == n, value]; if (length(x)) x[1] else default }
pseudocount <- as.numeric(param("pseudocount", "0.1"))
duplicate_gse <- param("duplicate_gse", "GSE253759")

safe_mkdir(file.path(PROJ, "results", "time_response"))
safe_mkdir(file.path(PROJ, "figures", "time_response"))

time_mids <- c("GSE215079", "GSE245014_Rab7", "GSE245014_Vps35", "GSE322548_TPM")
time_manifest <- data.frame(
  matrix_id = time_mids,
  gse = c("GSE215079", "GSE245014", "GSE245014", "GSE322548"),
  biological_design = c("AoAge2 deletion", "Rab7 deletion", "Vps35 deletion", "Hog1 knockdown"),
  independent_gse = TRUE,
  shared_wt_with_other_matrix = c(FALSE, TRUE, TRUE, FALSE),
  time_points = c("0,6,12,24", "0,12,36", "0,12,36", "0,12,24"),
  stringsAsFactors = FALSE
)
write_tsv(time_manifest, file.path(PROJ, "results", "time_response", "time_series_manifest.tsv"))

## Candidate trajectories are calculated inside each matrix in log2(TPM/FPKM
## + pseudocount) space. No values from different GEO records are concatenated.
x <- expr[matrix_id %in% time_mids & unit != "count" & is.finite(value)]
x[, value_log2 := log2(value + pseudocount)]
traj <- x[, .(n = sum(is.finite(value_log2)), mean_log2 = if (any(is.finite(value_log2))) mean(value_log2, na.rm = TRUE) else NA_real_,
              sd_log2 = if (sum(is.finite(value_log2)) >= 2L) sd(value_log2, na.rm = TRUE) else NA_real_,
              mean_value = if (any(is.finite(value))) mean(value, na.rm = TRUE) else NA_real_),
           by = .(candidate_id, matrix_id, gse, genotype, time_h, time_label, treatment)]
base <- traj[is.finite(time_h) & time_h == 0, .(candidate_id, matrix_id, genotype, baseline_mean_log2 = mean_log2, baseline_n = n)]
traj <- merge(traj, base, by = c("candidate_id", "matrix_id", "genotype"), all.x = TRUE)
traj[, induction_log2fc := mean_log2 - baseline_mean_log2]
traj[, induction_se := sd_log2 / sqrt(pmax(n, 1))]
traj[, induction_ci_low := induction_log2fc - 1.96 * induction_se]
traj[, induction_ci_high := induction_log2fc + 1.96 * induction_se]
traj[, direction := fifelse(induction_log2fc > 0, "up", fifelse(induction_log2fc < 0, "down", "flat"))]
setorder(traj, matrix_id, candidate_id, genotype, time_h)
write_tsv(traj, file.path(PROJ, "results", "time_response", "candidate_time_trajectories.tsv"))

auc_trap <- function(t, y) {
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]
  if (length(t) < 2L) return(NA_real_)
  o <- order(t); t <- t[o]; y <- y[o]
  sum(diff(t) * (head(y, -1L) + tail(y, -1L)) / 2)
}

classify_track <- function(d) {
  d <- d[is.finite(time_h) & time_h > 0 & is.finite(induction_log2fc)]
  if (!nrow(d)) return(data.table(n_time = 0L, max_abs_log2fc = NA_real_, peak_time_h = NA_real_, endpoint_log2fc = NA_real_,
                                  auc_signed = NA_real_, auc_abs = NA_real_, direction_consistency = NA_real_, trajectory_class = "not_detected"))
  setorder(d, time_h)
  y <- d$induction_log2fc; tt <- d$time_h; mx <- max(abs(y), na.rm = TRUE); peak <- tt[which.max(abs(y))][1]; endpoint <- y[which.max(tt)]
  pos <- sum(y > 0); neg <- sum(y < 0); cons <- max(pos, neg) / max(1, pos + neg)
  early_peak <- peak < max(tt) && abs(endpoint) <= 0.5 * mx && (sign(endpoint) == sign(y[which.max(abs(y))]) || abs(endpoint) < 0.2)
  late_peak <- peak == max(tt) && (length(y) == 1L || max(abs(y[tt < peak]), na.rm = TRUE) <= 0.5 * mx)
  persistent <- cons >= 0.75 && abs(endpoint) >= 0.7 * mx
  cls <- if (mx < 0.5 || any(d$n < 2L)) "nonresponse_or_low" else if (early_peak) "early" else if (late_peak) "late" else if (persistent) "persistent" else "complex"
  data.table(n_time = nrow(d), max_abs_log2fc = mx, peak_time_h = peak, endpoint_log2fc = endpoint,
             auc_signed = auc_trap(tt, y), auc_abs = auc_trap(tt, abs(y)), direction_consistency = cons, trajectory_class = cls)
}

## WT trajectories define the time-response phenotype; mutant trajectories and
## interaction tracks are retained in parallel for pathway-dependence calls.
wt_class <- traj[genotype == "WT", classify_track(.SD), by = .(candidate_id, matrix_id, gse)]
setnames(wt_class, c("trajectory_class", "max_abs_log2fc", "peak_time_h"), c("wt_trajectory_class", "wt_max_abs_log2fc", "wt_peak_time_h"), skip_absent = TRUE)
mut_class <- traj[genotype != "WT", classify_track(.SD), by = .(candidate_id, matrix_id, gse, genotype)]
setnames(mut_class, c("trajectory_class", "max_abs_log2fc", "peak_time_h"), c("mutant_trajectory_class", "mutant_max_abs_log2fc", "mutant_peak_time_h"), skip_absent = TRUE)
write_tsv(wt_class, file.path(PROJ, "results", "time_response", "wt_time_classification_by_matrix.tsv"))
write_tsv(mut_class, file.path(PROJ, "results", "time_response", "mutant_time_classification_by_matrix.tsv"))

interaction <- models[grepl("genotype_time_interaction", contrast_type) & matrix_id %in% time_mids & is.finite(log2fc),
                      .(candidate_id, matrix_id, gse, mutant, time_h, interaction_log2fc = log2fc, interaction_p_adj = p_adj, interaction_supported = supported_strong)]
interaction[, induction_log2fc := interaction_log2fc]
int_class <- interaction[, {
  classify_track(.SD[, .(time_h, induction_log2fc, n = 3L)])
}, by = .(candidate_id, matrix_id, gse, mutant)]
setnames(int_class, c("trajectory_class", "max_abs_log2fc", "peak_time_h"), c("interaction_trajectory_class", "interaction_max_abs_log2fc", "interaction_peak_time_h"), skip_absent = TRUE)
write_tsv(interaction, file.path(PROJ, "results", "time_response", "candidate_interaction_trajectories.tsv"))
write_tsv(int_class, file.path(PROJ, "results", "time_response", "interaction_time_classification_by_matrix.tsv"))

## Consensus is counted by distinct GEO accession, so Rab7 and Vps35 do not
## create two independent WT trajectory votes for GSE245014.
cons_gse <- wt_class[, .SD[which.max(wt_max_abs_log2fc)], by = .(candidate_id, gse)]
cons <- cons_gse[, {
  cls <- wt_trajectory_class
  tab <- sort(table(cls), decreasing = TRUE)
  top <- names(tab)[1]
  data.table(n_independent_gse = .N, classes = paste(paste(names(tab), as.integer(tab), sep = ":"), collapse = ";"),
             consensus_class = if (tab[1] >= 2L) top else "mixed_or_single_series",
             consensus_peak_time_h = if (top %in% cls) as.numeric(median(as.numeric(wt_peak_time_h[cls == top]), na.rm = TRUE)) else NA_real_,
             consensus_max_abs_log2fc = if (top %in% cls) as.numeric(max(as.numeric(wt_max_abs_log2fc[cls == top]), na.rm = TRUE)) else NA_real_)
}, by = candidate_id]
cons <- merge(data.table(candidate_id = candidates$sequence_id), cons, by = "candidate_id", all.x = TRUE)
write_tsv(cons, file.path(PROJ, "results", "time_response", "wt_time_consensus.tsv"))

## Compact trajectory heatmaps for rapid review.
plot_traj <- traj[genotype == "WT" & is.finite(induction_log2fc), .(log2fc = mean(induction_log2fc)), by = .(candidate_id, matrix_id, time_h)]
for (mid in time_mids) {
  h <- plot_traj[matrix_id == mid]
  if (!nrow(h)) next
  hm <- dcast(h, candidate_id ~ time_h, value.var = "log2fc")
  mat <- as.matrix(hm[, -1, with = FALSE]); rownames(mat) <- hm$candidate_id
  png(file.path(PROJ, "figures", "time_response", paste0(mid, "_WT_trajectory_heatmap.png")), width = 1800, height = 2600, res = 180)
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, na_col = "grey90", main = paste(mid, "WT induction trajectory"))
  dev.off()
}
top <- cons[order(-consensus_max_abs_log2fc)][1:min(12L, .N)]$candidate_id
if (length(top)) {
  p <- traj[genotype == "WT" & candidate_id %in% top & is.finite(induction_log2fc)]
  png(file.path(PROJ, "figures", "time_response", "top_WT_trajectories.png"), width = 1800, height = 1200, res = 180)
  print(ggplot(p, aes(time_h, induction_log2fc, group = candidate_id, color = candidate_id)) + geom_hline(yintercept = 0, linetype = 2) + geom_line() + geom_point() + facet_wrap(~ matrix_id, scales = "free_x") + theme_bw(base_size = 9) + labs(x = "hours", y = "WT induction log2FC", color = "candidate"))
  dev.off()
}

out_files <- c(list.files(file.path(PROJ, "results", "time_response"), full.names = TRUE), list.files(file.path(PROJ, "figures", "time_response"), full.names = TRUE))
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE),
          file.path(PROJ, "checksums", "step08_time_response_sha256.tsv"))
cat("time_series_matrices:", length(time_mids), "\n")
cat("trajectory_rows:", nrow(traj), "\n")
cat("wt_classified_rows:", nrow(wt_class), "\n")
cat("consensus_candidates:", nrow(cons), "\n")
cat("consensus_class_counts:\n"); print(cons[, .N, by = consensus_class])
finish_log()
