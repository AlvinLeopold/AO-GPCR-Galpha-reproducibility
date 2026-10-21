#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(DESeq2)
  library(ggplot2)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "07_series_models"
input_paths <- c(
  file.path(PROJ, "results", "expression", "parsed_matrices.rds"),
  file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv"),
  file.path(PROJ, "results", "design", "contrast_manifest.tsv"),
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "config", "deep_analysis_parameters.tsv")
)
start_log(step, input_paths)

parsed <- readRDS(input_paths[1])
samples <- fread(input_paths[2], sep = "\t", encoding = "UTF-8")
contrasts <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
candidates <- fread(input_paths[4], sep = "\t", encoding = "UTF-8")
params <- fread(input_paths[5], sep = "\t", encoding = "UTF-8")
param <- function(n, default) { x <- params[parameter == n, value]; if (length(x)) x[1] else default }
pseudocount <- as.numeric(param("pseudocount", "0.1"))
fdr_cut <- as.numeric(param("fdr_cutoff", "0.05"))
lfc_cut <- as.numeric(param("abs_log2fc_cutoff", "1"))
count_min_total <- as.numeric(param("count_min_total", "10"))
count_min_samples <- as.numeric(param("count_min_samples", "2"))
duplicate_gse <- param("duplicate_gse", "GSE253759")
order_gse <- param("order_based_gse", "GSE303416")

safe_mkdir(file.path(PROJ, "results", "models"))
safe_mkdir(file.path(PROJ, "results", "models", "full_count_results"))
safe_mkdir(file.path(PROJ, "figures", "model_qc"))

group_key <- function(genotype, time_h) paste(genotype, ifelse(is.finite(time_h), paste0(time_h, "h"), "Day5"), sep = "__")
coef_name <- function(x) paste0("g_", make.names(x))
contrast_weights <- function(ct, group_levels) {
  w <- setNames(rep(0, length(group_levels)), group_levels)
  mut <- ct$mutant; tt <- ct$time_h
  set <- function(key, value) { if (!key %in% names(w)) stop("Missing contrast group: ", key); w[key] <<- w[key] + value }
  if (ct$contrast_type == "wt_induction") { set(group_key("WT", tt), 1); set(group_key("WT", 0), -1)
  } else if (ct$contrast_type == "mutant_induction") { set(group_key(mut, tt), 1); set(group_key(mut, 0), -1)
  } else if (ct$contrast_type == "genotype_effect") {
    set(group_key(mut, tt), 1); set(group_key("WT", tt), -1)
  } else if (ct$contrast_type == "genotype_time_interaction") {
    set(group_key(mut, tt), 1); set(group_key(mut, 0), -1); set(group_key("WT", tt), -1); set(group_key("WT", 0), 1)
  } else stop("Unsupported contrast type: ", ct$contrast_type)
  w
}

candidate_map <- function(rec, gene_names) {
  rm <- as.data.table(rec$rowmap)
  rm[, model_gene_id := rownames(rec$values)[row_index]]
  rm <- rm[!is.na(candidate_id)]
  rm[match(gene_names, model_gene_id), .(candidate_id, raw_gene_id, raw_gene_name, mapping_method)]
}

result_rows <- list(); status_rows <- list(); not_tested_rows <- list(); pca_rows <- list(); full_count_files <- character()
for (mid in names(parsed)) {
  cat("modeling", mid, "\n")
  rec <- parsed[[mid]]; sm <- copy(samples[matrix_id == mid]); cm <- copy(contrasts[matrix_id == mid])
  mat <- rec$values
  if (!identical(as.character(rec$columns), as.character(sm$matrix_column))) stop("Sample order mismatch: ", mid)
  group <- factor(group_key(sm$genotype, sm$time_h), levels = unique(group_key(sm$genotype, sm$time_h)))
  design <- model.matrix(~0 + group); colnames(design) <- coef_name(levels(group))
  rownames(design) <- sm$sample_key
  keep <- rep(TRUE, nrow(mat)); method <- if (rec$unit == "count") "DESeq2" else "limma_log2"
  model_role <- if (rec$unit == "count") "confirmatory_raw_counts" else "discovery_processed_expression"
  status <- "success"; message <- ""; fit_obj <- NULL; model_mat <- NULL
  tryCatch({
    if (rec$unit == "count") {
      keep <- rowSums(mat, na.rm = TRUE) >= count_min_total & rowSums(mat > 0, na.rm = TRUE) >= count_min_samples &
        apply(mat, 1, function(z) all(is.finite(z)))
      count_mat <- round(mat[keep, , drop = FALSE])
      ## parsed matrices use source-file column names; DESeq2 requires exact
      ## identity with colData row names, so replace them with stable sample keys.
      colnames(count_mat) <- sm$sample_key
      dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = data.frame(row.names = sm$sample_key, group = group), design = ~0 + group)
      dds <- DESeq(dds, quiet = TRUE)
      fit_obj <- dds
      model_mat <- assay(vst(dds, blind = FALSE))
    } else {
      ## A few supplementary tables contain nonnumeric placeholders in genes
      ## outside the candidate set. Exclude such genes from limma rather than
      ## silently converting a missing value into zero.
      keep <- apply(mat, 1, function(z) all(is.finite(z)))
      model_mat <- log2(mat[keep, , drop = FALSE] + pseudocount)
      fit_obj <- eBayes(lmFit(model_mat, design), trend = TRUE, robust = TRUE)
    }
  }, error = function(e) { status <<- "failed"; message <<- conditionMessage(e) })

  if (status == "success") {
    pc <- prcomp(t(model_mat), center = TRUE, scale. = FALSE)
    var_exp <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    pca_rows[[length(pca_rows) + 1L]] <- data.frame(matrix_id = mid, gse = rec$gse, unit = rec$unit, sample_key = sm$sample_key, gsm = sm$gsm,
      genotype = sm$genotype, time_h = sm$time_h, PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC1_percent = var_exp[1], PC2_percent = var_exp[2], stringsAsFactors = FALSE)
    pp <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], genotype = sm$genotype, time_label = sm$time_label)
    png(file.path(PROJ, "figures", "model_qc", paste0(mid, "_PCA.png")), width = 1400, height = 1050, res = 160)
    print(ggplot(pp, aes(PC1, PC2, color = genotype, shape = time_label)) + geom_point(size = 3) + theme_bw(base_size = 10) + labs(title = paste(mid, rec$unit), shape = "time"))
    dev.off()
  }

  matrix_results <- list()
  if (status == "success") for (ii in seq_len(nrow(cm))) {
    ct <- cm[ii]
    one <- tryCatch({
      w <- contrast_weights(ct, levels(group))
      if (rec$unit == "count") {
        rn <- resultsNames(fit_obj); expected <- coef_name(levels(group))
        ## DESeq2 retains the factor name (group...) in coefficient names,
        ## whereas the contrast helper uses stable g_ names.
        actual <- sub("^g_", "group", expected)
        if (!all(actual %in% rn)) stop("DESeq2 coefficient mismatch: ", paste(setdiff(actual, rn), collapse = ","))
        cv <- setNames(rep(0, length(rn)), rn); cv[actual] <- w
        rr <- as.data.frame(results(fit_obj, contrast = cv, alpha = fdr_cut))
        rr$model_gene_id <- rownames(rr)
        setDT(rr); setnames(rr, c("log2FoldChange", "lfcSE", "pvalue", "padj"), c("log2fc", "se", "p_value", "p_adj"))
        rr[, `:=`(ci_low = log2fc - 1.96 * se, ci_high = log2fc + 1.96 * se)]
      } else {
        cv <- matrix(w, ncol = 1, dimnames = list(coef_name(levels(group)), "contrast"))
        cf <- contrasts.fit(fit_obj, cv); cf <- eBayes(cf, trend = TRUE, robust = TRUE)
        rr <- as.data.table(topTable(cf, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "model_gene_id")
        setnames(rr, c("logFC", "P.Value", "adj.P.Val"), c("log2fc", "p_value", "p_adj"))
        ## limma's moderated t statistic gives an approximate standard error
        ## for the log2 effect; retain it for confidence intervals.
        if (!"t" %in% names(rr)) rr[, t := NA_real_]
        rr[, se := ifelse(is.finite(get("t")) & get("t") != 0, abs(log2fc / get("t")), NA_real_)]
        rr[!is.finite(se), se := NA_real_]
        rr[, `:=`(ci_low = log2fc - 1.96 * se, ci_high = log2fc + 1.96 * se)]
      }
      rr[, `:=`(matrix_id = mid, gse = rec$gse, unit = rec$unit, method = method, model_role = model_role,
                 contrast_id = ct$contrast_id, contrast_type = ct$contrast_type, mutant = ct$mutant, time_h = ct$time_h,
                 mapping_confidence = paste(unique(sm$mapping_confidence), collapse = ";"),
                 independent_evidence = rec$gse != duplicate_gse, provisional_mapping = rec$gse == order_gse)]
      rr[, supported_fdr := is.finite(p_adj) & p_adj < fdr_cut]
      rr[, supported_strong := supported_fdr & is.finite(log2fc) & abs(log2fc) >= lfc_cut]
      rr
    }, error = function(e) {
      not_tested_rows[[length(not_tested_rows) + 1L]] <<- data.frame(matrix_id = mid, gse = rec$gse, unit = rec$unit, contrast_id = ct$contrast_id,
        contrast_type = ct$contrast_type, mutant = ct$mutant, time_h = ct$time_h, reason = conditionMessage(e), stringsAsFactors = FALSE)
      NULL
    })
    if (!is.null(one)) matrix_results[[length(matrix_results) + 1L]] <- one
  }
  all_mid <- rbindlist(matrix_results, fill = TRUE)
  if (nrow(all_mid)) {
    if (rec$unit == "count") {
      fpath <- file.path(PROJ, "results", "models", "full_count_results", paste0(mid, "_all_genes.tsv.gz"))
      fwrite(all_mid, fpath, sep = "\t", quote = FALSE, na = "NA"); full_count_files <- c(full_count_files, fpath)
    }
    cmapping <- candidate_map(rec, all_mid$model_gene_id)
    cand <- cbind(all_mid, cmapping)
    cand <- cand[!is.na(candidate_id)]
    result_rows[[length(result_rows) + 1L]] <- cand
  }
  status_rows[[length(status_rows) + 1L]] <- data.frame(matrix_id = mid, gse = rec$gse, unit = rec$unit, method = method, model_role = model_role,
    status = status, message = message, n_samples = nrow(sm), n_genes = nrow(mat), n_genes_tested = sum(keep), n_contrasts_planned = nrow(cm),
    n_contrasts_completed = length(matrix_results), mapping_confidence = paste(unique(sm$mapping_confidence), collapse = ";"), stringsAsFactors = FALSE)
}

res <- rbindlist(result_rows, fill = TRUE)
setcolorder(res, c("candidate_id", "matrix_id", "gse", "unit", "method", "model_role", "contrast_id", "contrast_type", "mutant", "time_h",
                   "log2fc", "se", "ci_low", "ci_high", "p_value", "p_adj", "supported_fdr", "supported_strong", "independent_evidence",
                   "provisional_mapping", "mapping_confidence", "model_gene_id", "raw_gene_id", "raw_gene_name", "mapping_method"))
write_tsv(res, file.path(PROJ, "results", "models", "candidate_series_model_results.tsv"))
write_tsv(res[supported_fdr == TRUE], file.path(PROJ, "results", "models", "candidate_series_model_significant.tsv"))
write_tsv(rbindlist(status_rows, fill = TRUE), file.path(PROJ, "results", "models", "model_status.tsv"))
not_tested <- if (length(not_tested_rows)) rbindlist(not_tested_rows, fill = TRUE) else data.table(
  matrix_id = character(), gse = character(), unit = character(), contrast_id = character(), contrast_type = character(),
  mutant = character(), time_h = numeric(), reason = character())
write_tsv(not_tested, file.path(PROJ, "results", "models", "series_models.not_tested.tsv"))
write_tsv(rbindlist(pca_rows, fill = TRUE), file.path(PROJ, "results", "models", "sample_pca_coordinates.tsv"))

## Direct concordance between paired processed and raw-count result matrices.
pair_map <- data.table(processed = c("GSE301423_TPM", "GSE301452_FPKM", "GSE303416_FPKM", "GSE306156_TPM", "GSE322548_TPM"),
                       counts = c("GSE301423_counts", "GSE301452_counts", "GSE303416_counts", "GSE306156_counts", "GSE322548_counts"))
conc_rows <- list()
for (i in seq_len(nrow(pair_map))) {
  a <- res[matrix_id == pair_map$processed[i], .(candidate_id, contrast_type, mutant, time_h, lfc_processed = log2fc, fdr_processed = p_adj)]
  b <- res[matrix_id == pair_map$counts[i], .(candidate_id, contrast_type, mutant, time_h, lfc_counts = log2fc, fdr_counts = p_adj)]
  z <- merge(a, b, by = c("candidate_id", "contrast_type", "mutant", "time_h"))
  z[, `:=`(processed_matrix = pair_map$processed[i], count_matrix = pair_map$counts[i], direction_agree = sign(lfc_processed) == sign(lfc_counts))]
  conc_rows[[length(conc_rows) + 1L]] <- z
}
concordance <- rbindlist(conc_rows, fill = TRUE)
write_tsv(concordance, file.path(PROJ, "results", "models", "processed_count_concordance.tsv"))

out_files <- c(list.files(file.path(PROJ, "results", "models"), full.names = TRUE, recursive = TRUE),
               list.files(file.path(PROJ, "figures", "model_qc"), full.names = TRUE))
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE),
          file.path(PROJ, "checksums", "step07_models_sha256.tsv"))
cat("model_matrices:", nrow(rbindlist(status_rows)), "\n")
cat("successful_matrices:", sum(rbindlist(status_rows)$status == "success"), "\n")
cat("candidate_result_rows:", nrow(res), "\n")
cat("candidate_fdr_rows:", sum(res$supported_fdr, na.rm = TRUE), "\n")
cat("candidate_strong_rows:", sum(res$supported_strong, na.rm = TRUE), "\n")
cat("not_tested_rows:", nrow(not_tested), "\n")
cat("paired_direction_agreement:", mean(concordance$direction_agree, na.rm = TRUE), "\n")
finish_log()
