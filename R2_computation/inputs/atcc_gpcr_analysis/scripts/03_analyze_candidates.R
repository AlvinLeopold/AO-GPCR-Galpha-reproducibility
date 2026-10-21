#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "03_analyze_candidates"
long_path <- file.path(PROJ, "results", "expression", "candidate_expression_long.tsv")
sample_path <- file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv")
parsed_path <- file.path(PROJ, "results", "expression", "parsed_matrices.rds")
candidates_path <- file.path(PROJ, "metadata", "candidates.tsv")
start_log(step, c(long_path, sample_path, parsed_path, candidates_path))

long <- fread(long_path, sep = "\t", encoding = "UTF-8")
samples <- fread(sample_path, sep = "\t", encoding = "UTF-8")
candidates <- fread(candidates_path, sep = "\t", encoding = "UTF-8")
parsed <- readRDS(parsed_path)

safe_mean <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_median <- function(x) if (all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
safe_cv <- function(x) { m <- safe_mean(x); s <- safe_sd(x); if (!is.finite(m) || m == 0 || !is.finite(s)) NA_real_ else abs(s / m) }
safe_ttest <- function(a, b, paired = FALSE) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  if (paired && length(a) == length(b)) return(tryCatch(t.test(a, b, paired = TRUE)$p.value, error = function(e) NA_real_))
  tryCatch(t.test(a, b, var.equal = FALSE)$p.value, error = function(e) NA_real_)
}

## Sample-level QC uses every gene in each parsed matrix, not only the 93 GPCRs.
qc_list <- list()
rep_qc_list <- list()
for (mid in names(parsed)) {
  rec <- parsed[[mid]]
  sm <- samples[matrix_id == mid]
  mat <- rec$values
  if (ncol(mat) != nrow(sm)) stop("Sample order mismatch in ", mid)
  sample_key <- sm$sample_key
  total <- colSums(mat, na.rm = TRUE)
  detected <- colSums(mat > 0, na.rm = TRUE)
  cand_sub <- long[matrix_id == mid]
  cand_det <- cand_sub[, .(candidate_detected = sum(detected_in_matrix & is.finite(value) & value > 0, na.rm = TRUE)), by = sample_key]
  sm2 <- merge(sm, cand_det, by = "sample_key", all.x = TRUE)
  sm2[is.na(candidate_detected), candidate_detected := 0L]
  lg <- log2(pmax(mat, 0) + 0.1)
  cor_med <- NA_real_
  if (ncol(lg) >= 2L) {
    cm <- suppressWarnings(cor(lg, use = "pairwise.complete.obs", method = "pearson"))
    cor_med <- safe_median(cm[upper.tri(cm)])
  }
  for (j in seq_len(nrow(sm2))) {
    qc_list[[length(qc_list) + 1L]] <- data.frame(gse = sm2$gse[j], matrix_id = mid, gsm = sm2$gsm[j],
      sample_key = sm2$sample_key[j], unit = rec$unit, genotype = sm2$genotype[j], time_h = sm2$time_h[j],
      replicate = sm2$replicate[j], total_signal_all_genes = total[j], detected_genes = detected[j],
      candidate_detected = sm2$candidate_detected[j], median_sample_correlation = cor_med,
      source_file = rec$file, stringsAsFactors = FALSE)
  }
  sub <- long[matrix_id == mid & unit != "count" & is.finite(value)]
  if (nrow(sub)) {
    z <- sub[, .(n = .N, mean_value = safe_mean(value), median_value = safe_median(value),
                 sd_value = safe_sd(value), cv = safe_cv(value), detection_rate = mean(value > 0)),
             by = .(candidate_id, gse, matrix_id, unit, genotype, time_h, time_label, treatment)]
    rep_qc_list[[length(rep_qc_list) + 1L]] <- z
  }
}
sample_qc <- rbindlist(qc_list, fill = TRUE)
rep_qc <- rbindlist(rep_qc_list, fill = TRUE)
write_tsv(sample_qc, file.path(PROJ, "results", "qc", "sample_qc.tsv"))
write_tsv(rep_qc, file.path(PROJ, "results", "qc", "candidate_replicate_qc.tsv"))

## Pairwise effects on normalized processed expression. A pseudo-count is used
## only for log-scale effect estimation and is recorded in the output.
pseudocount <- 0.1
norm_long <- long[unit != "count"]
effect_list <- list()
for (mid in unique(norm_long$matrix_id)) {
  subm <- norm_long[matrix_id == mid & is.finite(value)]
  if (!nrow(subm)) next
  groups <- unique(subm[, .(genotype, time_h, time_label)])
  groups <- groups[order(genotype, time_h)]
  genotypes <- unique(groups$genotype)
  for (cid in candidates$sequence_id) {
    sg <- subm[candidate_id == cid]
    if (!nrow(sg)) next
    # Within-genotype induction effects, with time zero as baseline.
    for (g in genotypes) {
      gg <- sg[genotype == g]
      if (!nrow(gg)) next
      times <- sort(unique(gg$time_h[is.finite(gg$time_h)]))
      if (!any(times == 0)) next
      for (tt in times[times > 0]) {
        a <- gg[time_h == 0, log2(value + pseudocount)]
        b <- gg[time_h == tt, log2(value + pseudocount)]
        if (length(a) < 2L || length(b) < 2L) next
        effect_list[[length(effect_list) + 1L]] <- data.frame(candidate_id = cid, matrix_id = mid, gse = unique(sg$gse), unit = unique(sg$unit),
          contrast_type = if (g == "WT") "wt_induction" else "mutant_induction", genotype = g, time_h = tt,
          group_a = paste0(g, "__0h"), group_b = paste0(g, "__", tt, "h"), n_a = length(a), n_b = length(b),
          mean_a = mean(a), mean_b = mean(b), log2fc = mean(b) - mean(a), pvalue = safe_ttest(a, b),
          pseudocount = pseudocount, stringsAsFactors = FALSE)
      }
    }
    # Mutant-vs-WT effects at each time.
    mutants <- setdiff(genotypes, "WT")
    for (mut in mutants) {
      mt <- sg[genotype == mut]
      wt <- sg[genotype == "WT"]
      if (!nrow(mt) || !nrow(wt)) next
      times <- intersect(unique(mt$time_h[is.finite(mt$time_h)]), unique(wt$time_h[is.finite(wt$time_h)]))
      for (tt in sort(times)) {
        a <- wt[time_h == tt, log2(value + pseudocount)]
        b <- mt[time_h == tt, log2(value + pseudocount)]
        if (length(a) < 2L || length(b) < 2L) next
        effect_list[[length(effect_list) + 1L]] <- data.frame(candidate_id = cid, matrix_id = mid, gse = unique(sg$gse), unit = unique(sg$unit),
          contrast_type = "genotype_effect", genotype = mut, time_h = tt, group_a = paste0("WT__", tt, "h"),
          group_b = paste0(mut, "__", tt, "h"), n_a = length(a), n_b = length(b), mean_a = mean(a), mean_b = mean(b),
          log2fc = mean(b) - mean(a), pvalue = safe_ttest(a, b), pseudocount = pseudocount, stringsAsFactors = FALSE)
      }
    }
    # Interaction proxy: paired changes in mutant and WT when replicate counts match.
    for (mut in mutants) {
      mt <- sg[genotype == mut]; wt <- sg[genotype == "WT"]
      if (!nrow(mt) || !nrow(wt)) next
      times <- intersect(unique(mt$time_h[is.finite(mt$time_h)]), unique(wt$time_h[is.finite(wt$time_h)]))
      for (tt in sort(times[times > 0])) {
        m0 <- mt[time_h == 0][order(replicate), log2(value + pseudocount)]
        mtv <- mt[time_h == tt][order(replicate), log2(value + pseudocount)]
        w0 <- wt[time_h == 0][order(replicate), log2(value + pseudocount)]
        wtv <- wt[time_h == tt][order(replicate), log2(value + pseudocount)]
        n <- min(length(m0), length(mtv), length(w0), length(wtv))
        if (n < 2L) next
        dm <- mtv[seq_len(n)] - m0[seq_len(n)]
        dw <- wtv[seq_len(n)] - w0[seq_len(n)]
        effect_list[[length(effect_list) + 1L]] <- data.frame(candidate_id = cid, matrix_id = mid, gse = unique(sg$gse), unit = unique(sg$unit),
          contrast_type = "genotype_time_interaction", genotype = mut, time_h = tt, group_a = paste0("WT_delta_", tt, "h"),
          group_b = paste0(mut, "_delta_", tt, "h"), n_a = n, n_b = n, mean_a = mean(dw), mean_b = mean(dm),
          log2fc = mean(dm) - mean(dw), pvalue = safe_ttest(dw, dm, paired = TRUE), pseudocount = pseudocount, stringsAsFactors = FALSE)
      }
    }
    # Day-5/one-state experiments still contribute genotype effects.
  }
}
effects <- if (length(effect_list)) rbindlist(effect_list, fill = TRUE) else data.table()
if (nrow(effects)) {
  effects[, p_adj := p.adjust(pvalue, method = "BH"), by = .(matrix_id, contrast_type, time_h)]
  effects[, supported := is.finite(p_adj) & p_adj < 0.05 & abs(log2fc) >= 1]
}
write_tsv(effects, file.path(PROJ, "results", "effects", "candidate_effects.tsv"))

## Baseline summaries are kept within each matrix/unit and never averaged across
## TPM and FPKM. The integration step uses counts of independent matrix records.
baseline <- norm_long[baseline == TRUE & genotype == "WT" & is.finite(value),
  .(n = .N, mean_value = safe_mean(value), median_value = safe_median(value), mean_log2 = safe_mean(log2(value + pseudocount)),
    detection_rate = mean(value > 0), replicate_cv = safe_cv(value)), by = .(candidate_id, gse, matrix_id, unit, time_h)]
write_tsv(baseline, file.path(PROJ, "results", "effects", "baseline_expression.tsv"))

## DESeq2 on the five count-bearing matrices. Pairwise group contrasts are
## generated for all groups with >=2 biological replicates; interaction contrasts
## are also emitted as a list contrast where the design permits it.
count_ids <- unique(long[unit == "count", matrix_id])
deseq_list <- list(); deseq_status <- list()
for (mid in count_ids) {
  rec <- parsed[[mid]]
  sm <- samples[matrix_id == mid]
  mat <- rec$values
  mat[!is.finite(mat)] <- 0
  mat <- round(pmax(mat, 0))
  rownames(mat) <- make.unique(rownames(mat))
  colnames(mat) <- sm$sample_key
  group <- paste(sm$genotype, ifelse(is.finite(sm$time_h), sm$time_h, "NA"), sep = "__")
  keep <- rowSums(mat) >= 10
  status <- "success"
  err <- ""
  if (sum(keep) < 10L || length(unique(group)) < 2L || any(table(group) < 2L)) {
    status <- "not_tested_insufficient_groups"
    err <- "Need at least two groups and two replicates per group after count filtering"
  } else {
    dds <- tryCatch({
      d <- DESeqDataSetFromMatrix(countData = mat[keep, , drop = FALSE],
                                  colData = data.frame(group = factor(group), row.names = colnames(mat)), design = ~group)
      d <- DESeq(d, quiet = TRUE)
      d
    }, error = function(e) { status <<- "failed"; err <<- conditionMessage(e); NULL })
    if (!is.null(dds)) {
      lev <- levels(colData(dds)$group)
      contrast_rows <- list()
      # All ordered pairwise group comparisons.
      for (a in lev) for (b in lev) if (a != b && match(a, lev) < match(b, lev)) {
        rr <- tryCatch(as.data.frame(results(dds, contrast = c("group", b, a))), error = function(e) NULL)
        if (is.null(rr)) next
        rr$raw_gene_id <- rownames(rr)
        rr$matrix_id <- mid; rr$gse <- rec$gse; rr$unit <- rec$unit
        rr$contrast_type <- "DESeq2_pairwise"; rr$group_a <- a; rr$group_b <- b
        contrast_rows[[length(contrast_rows) + 1L]] <- rr
      }
      # Difference-in-differences contrast for groups sharing WT/mutant and time.
      muts <- setdiff(unique(sm$genotype), "WT")
      for (mut in muts) {
        ts <- sort(unique(sm$time_h[is.finite(sm$time_h)]))
        for (tt in ts[ts > 0]) {
          needed <- c(paste0("WT__0"), paste0("WT__", tt), paste0(mut, "__0"), paste0(mut, "__", tt))
          if (!all(needed %in% lev)) next
          rr <- tryCatch(as.data.frame(results(dds, contrast = list(c("group", paste0(mut, "__", tt), paste0("WT__", tt)),
                                                               c("group", paste0(mut, "__0"), paste0("WT__0"))))), error = function(e) NULL)
          if (is.null(rr)) next
          rr$raw_gene_id <- rownames(rr); rr$matrix_id <- mid; rr$gse <- rec$gse; rr$unit <- rec$unit
          rr$contrast_type <- "DESeq2_interaction"; rr$group_a <- paste0("WT_delta_", tt); rr$group_b <- paste0(mut, "_delta_", tt)
          contrast_rows[[length(contrast_rows) + 1L]] <- rr
        }
      }
      if (length(contrast_rows)) {
        rr <- rbindlist(contrast_rows, fill = TRUE)
        rowmap <- rec$rowmap
        rr$candidate_id <- NA_character_
        m <- match(rr$raw_gene_id, make.unique(as.character(rowmap$raw_gene_id)))
        rr$candidate_id <- rowmap$candidate_id[m]
        rr <- rr[!is.na(candidate_id)]
        if (nrow(rr)) deseq_list[[length(deseq_list) + 1L]] <- rr
      }
    }
  }
  deseq_status[[length(deseq_status) + 1L]] <- data.frame(matrix_id = mid, gse = rec$gse, unit = rec$unit,
    status = status, message = err, n_samples = nrow(sm), n_groups = length(unique(group)),
    n_genes = nrow(mat), n_genes_after_filter = sum(keep), stringsAsFactors = FALSE)
}
deseq <- if (length(deseq_list)) rbindlist(deseq_list, fill = TRUE) else data.table()
if (nrow(deseq)) {
  setcolorder(deseq, c("candidate_id", "matrix_id", "gse", "unit", "contrast_type", "group_a", "group_b", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "raw_gene_id"))
  deseq[, supported := is.finite(padj) & padj < 0.05 & abs(log2FoldChange) >= 1]
}
write_tsv(deseq, file.path(PROJ, "results", "deseq2", "deseq2_candidate_results.tsv"))
write_tsv(rbindlist(deseq_status, fill = TRUE), file.path(PROJ, "results", "deseq2", "deseq2_status.tsv"))

## Compact candidate-level score components for the integration step.
base_summary <- baseline[, .(baseline_matrix_count = uniqueN(matrix_id), baseline_detected_matrix_count = sum(detection_rate > 0),
                             baseline_mean_log2 = safe_mean(mean_log2), baseline_max_log2 = max(mean_log2, na.rm = TRUE)), by = candidate_id]
if (nrow(effects)) {
  eff_summary <- effects[, .(effect_records = .N, supported_effect_records = sum(supported, na.rm = TRUE),
                             induction_records = sum(contrast_type %in% c("wt_induction", "mutant_induction")),
                             induction_supported = sum(contrast_type %in% c("wt_induction", "mutant_induction") & supported, na.rm = TRUE),
                             genotype_supported = sum(contrast_type == "genotype_effect" & supported, na.rm = TRUE),
                             interaction_supported = sum(contrast_type == "genotype_time_interaction" & supported, na.rm = TRUE)), by = candidate_id]
} else eff_summary <- data.table()
if (nrow(deseq)) {
  de_summary <- deseq[, .(deseq_records = .N, deseq_supported = sum(supported, na.rm = TRUE),
                          deseq_interaction_supported = sum(contrast_type == "DESeq2_interaction" & supported, na.rm = TRUE)), by = candidate_id]
} else de_summary <- data.table()
score <- merge(candidates[, .(candidate_id = sequence_id, sequence_id, tm_count, gpcr_domain_count, gpcr_domains, candidate_order)], base_summary, by = "candidate_id", all.x = TRUE)
score <- merge(score, eff_summary, by = "candidate_id", all.x = TRUE)
score <- merge(score, de_summary, by = "candidate_id", all.x = TRUE)
for (j in setdiff(colnames(score), c("candidate_id", "sequence_id", "gpcr_domains"))) if (is.numeric(score[[j]])) score[[j]][!is.finite(score[[j]])] <- 0
score[is.na(score)] <- 0
write_tsv(score, file.path(PROJ, "results", "effects", "candidate_score_components.tsv"))

## Figures.
safe_mkdir(file.path(PROJ, "figures"))
png(file.path(PROJ, "figures", "sample_library_size.png"), width = 1800, height = 1200, res = 180)
print(ggplot(sample_qc, aes(x = matrix_id, y = log10(total_signal_all_genes + 1), fill = unit)) +
        geom_boxplot(outlier.size = 0.5) + coord_flip() + theme_bw(base_size = 8) +
        labs(x = NULL, y = "log10(total signal + 1)", title = "GEO expression-table sample signal"))
dev.off()

if (nrow(baseline)) {
  bwide <- dcast(baseline, candidate_id ~ matrix_id, value.var = "mean_log2")
  bm <- as.matrix(bwide[, -1, with = FALSE]); rownames(bm) <- bwide$candidate_id
  bm[!is.finite(bm)] <- NA_real_
  png(file.path(PROJ, "figures", "candidate_baseline_heatmap.png"), width = 1900, height = 2600, res = 180)
  pheatmap(bm, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
           main = "WT baseline expression (log2(value + 0.1))", na_col = "grey90")
  dev.off()
}

score_plot <- copy(score)
score_plot[, ranking_score := scale01(baseline_mean_log2) + scale01(induction_supported) + scale01(genotype_supported) + scale01(interaction_supported) + scale01(deseq_supported)]
top <- score_plot[order(-ranking_score)][1:min(20L, .N)]
png(file.path(PROJ, "figures", "candidate_priority_top20.png"), width = 1800, height = 1200, res = 180)
print(ggplot(top, aes(x = reorder(candidate_id, ranking_score), y = ranking_score, fill = tm_count == 7)) +
        geom_col() + coord_flip() + theme_bw(base_size = 9) + theme(legend.position = "none") +
        labs(x = NULL, y = "evidence score (unweighted components)", title = "Top 20 candidates by preliminary evidence score"))
dev.off()

write_tsv(score_plot[order(-ranking_score)], file.path(PROJ, "results", "effects", "candidate_score_components_ranked.tsv"))
write_tsv(data.frame(path = list.files(file.path(PROJ, "figures"), full.names = TRUE),
                     sha256 = vapply(list.files(file.path(PROJ, "figures"), full.names = TRUE), file_sha256, character(1))),
          file.path(PROJ, "checksums", "step03_outputs_sha256.tsv"))

cat("normalized_effect_records:", nrow(effects), "\n")
cat("baseline_records:", nrow(baseline), "\n")
cat("sample_qc_records:", nrow(sample_qc), "\n")
cat("DESeq2_candidate_records:", nrow(deseq), "\n")
print(rbindlist(deseq_status, fill = TRUE))
finish_log()
