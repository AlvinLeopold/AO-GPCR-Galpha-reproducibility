#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(limma)
  library(fgsea)
  library(ggplot2)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)

step <- "10_pathway_network"
input_paths <- c(
  file.path(PROJ, "results", "expression", "parsed_matrices.rds"),
  file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv"),
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "results", "models", "candidate_series_model_results.tsv"),
  file.path(PROJ, "config", "deep_analysis_parameters.tsv")
)
start_log(step, input_paths)
parsed <- readRDS(input_paths[1]); samples <- fread(input_paths[2], sep = "\t", encoding = "UTF-8"); candidates <- fread(input_paths[3], sep = "\t", encoding = "UTF-8")
models <- fread(input_paths[4], sep = "\t", encoding = "UTF-8"); params <- fread(input_paths[5], sep = "\t", encoding = "UTF-8")
param <- function(n, default) { x <- params[parameter == n, value]; if (length(x)) x[1] else default }
fdr_cut <- as.numeric(param("fdr_cutoff", "0.05")); top_genes <- as.integer(param("network_top_variable_genes", "5000")); min_samples <- as.integer(param("network_min_samples", "15"))
min_module <- as.integer(param("network_min_module_size", "30")); cor_fdr <- as.numeric(param("correlation_fdr_cutoff", "0.05"))

safe_mkdir(file.path(PROJ, "results", "network")); safe_mkdir(file.path(PROJ, "figures", "network")); safe_mkdir(file.path(PROJ, "results", "pathway"))

## The installed environment has fgsea and WGCNA but not GSVA. Record the
## substitution explicitly and use module eigengenes plus fgsea where gene
## sets are available; no silent GSVA claim is made.
write_tsv(data.frame(component = c("GSVA", "WGCNA", "fgsea"), available = c(requireNamespace("GSVA", quietly = TRUE), TRUE, TRUE),
                         method_used = c("not_run_GSVA_missing_dependency", "signed_network_module_eigengenes", "ranked_gene_set_enrichment"), stringsAsFactors = FALSE),
          file.path(PROJ, "results", "network", "pathway_method_status.tsv"))

candidate_set <- candidates$sequence_id
network_mids <- names(parsed)[vapply(parsed, function(z) z$unit == "count" && ncol(z$values) >= min_samples, logical(1))]
network_manifest <- data.table(matrix_id = network_mids, gse = vapply(parsed[network_mids], `[[`, character(1), "gse"), unit = "count",
                               n_samples = vapply(parsed[network_mids], function(z) ncol(z$values), integer(1)), eligible = TRUE,
                               note = "raw counts; WGCNA exploratory and series-specific")
write_tsv(network_manifest, file.path(PROJ, "results", "network", "network_manifest.tsv"))

## Predefined fungal response gene sets. These are deliberately compact and
## auditable; users can replace the IDs with an annotation-derived GMT later.
gene_sets <- list(
  cAMP_PKA = c("AOL_s00003g", "AOL_s00004g"),
  MAPK_osmotic = c("AOL_s00001g", "AOL_s00002g"),
  autophagy = c("AOL_s00005g", "AOL_s00006g"),
  mitochondrial_fission = c("AOL_s00007g", "AOL_s00008g"),
  hyphal_fusion = c("AOL_s00009g", "AOL_s00010g"),
  secondary_metabolism = c("AOL_s00011g", "AOL_s00012g")
)
## The placeholder-style IDs above are not used for claims unless they match a
## real model gene ID; this keeps pathway output explicit rather than inventing
## annotations for unverified fungal loci.

module_rows <- list(); module_trait_rows <- list(); module_gene_rows <- list(); corr_rows <- list(); de_rows <- list()
for (mid in network_mids) {
  rec <- parsed[[mid]]; sm <- samples[matrix_id == mid]; mat <- rec$values
  keep <- rowSums(mat, na.rm = TRUE) >= 10 & rowSums(mat > 0, na.rm = TRUE) >= 2 & apply(mat, 1, function(z) all(is.finite(z)))
  mat <- mat[keep, , drop = FALSE]; genes <- rownames(mat)
  if (nrow(mat) < min_module * 3L) next
  vars <- apply(log2(mat + 1), 1, var); ord <- order(vars, decreasing = TRUE)[seq_len(min(top_genes, length(vars)))]
  datExpr <- t(log2(mat[ord, , drop = FALSE] + 1)); colnames(datExpr) <- genes[ord]
  qc <- goodSamplesGenes(datExpr, verbose = FALSE)
  datExpr <- datExpr[qc$goodSamples, qc$goodGenes, drop = FALSE]; sm2 <- sm[qc$goodSamples]
  if (nrow(datExpr) < min_samples || ncol(datExpr) < min_module * 3L) next
  powers <- c(1:10, 12, 14, 16, 18, 20); sft <- suppressWarnings(pickSoftThreshold(datExpr, powerVector = powers, verbose = 0))
  chosen <- 6L
  if (!is.null(sft$fitIndices)) {
    r2 <- sft$fitIndices[, 2]; ix <- which(r2 >= 0.8); if (length(ix)) chosen <- sft$fitIndices$Power[ix[1]]
  }
  adj <- adjacency(datExpr, power = chosen, type = "signed"); tom <- TOMsimilarity(adj, TOMType = "signed"); diss <- 1 - tom
  tree <- hclust(as.dist(diss), method = "average"); mods <- cutreeDynamic(dendro = tree, distM = diss, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min_module, verbose = 0)
  colors <- labels2colors(mods); eig <- moduleEigengenes(datExpr, colors)$eigengenes; eig <- orderMEs(eig)
  for (j in seq_along(unique(colors))) {
    col <- sort(unique(colors))[j]; genes_mod <- colnames(datExpr)[colors == col]; me_name <- paste0("ME", col)
    module_gene_rows[[length(module_gene_rows) + 1L]] <- data.frame(matrix_id = mid, gse = rec$gse, module = col, module_eigengene = me_name, gene = genes_mod, gene_is_candidate = genes_mod %in% normalize_id(candidate_set), stringsAsFactors = FALSE)
    module_rows[[length(module_rows) + 1L]] <- data.frame(matrix_id = mid, gse = rec$gse, module = col, module_eigengene = me_name, module_size = length(genes_mod), soft_power = chosen, n_samples = nrow(datExpr), stringsAsFactors = FALSE)
  }
  ## Candidate-module correlations use Spearman and BH within each matrix.
  ## Map candidate IDs through the parsed rowmap rather than assuming that
  ## model rownames are the published AOL identifiers (many matrices use
  ## gene### or MSTRG identifiers).
  candidate_model_ids <- rec$rowmap$model_gene_id <- rownames(rec$values)[rec$rowmap$row_index]
  for (cid in candidate_set) {
    ccore <- normalize_id(cid); hit <- which(normalize_id(rec$rowmap$candidate_id) == ccore & rec$rowmap$model_gene_id %in% colnames(datExpr))
    if (!length(hit)) next
    candidate_col <- rec$rowmap$model_gene_id[hit[1]]
    if (!candidate_col %in% colnames(datExpr)) next
    for (me in colnames(eig)) {
      a <- suppressWarnings(cor.test(datExpr[, candidate_col], eig[, me], method = "spearman", exact = FALSE))
      corr_rows[[length(corr_rows) + 1L]] <- data.frame(candidate_id = cid, matrix_id = mid, gse = rec$gse, module_eigengene = me,
        rho = unname(a$estimate), p_value = a$p.value, stringsAsFactors = FALSE)
    }
  }
  ## Module-level phenotype correlations support interpretation of perturbation
  ## effects without merging samples across GEO studies.
  ph <- data.frame(induction = as.numeric(sm2$time_h > 0), mutant = as.numeric(sm2$genotype != "WT"), row.names = sm2$sample_key)
  for (me in colnames(eig)) {
    for (trait in colnames(ph)) {
      a <- suppressWarnings(cor.test(eig[, me], ph[[trait]], method = "spearman", exact = FALSE))
      module_trait_rows[[length(module_trait_rows) + 1L]] <- data.frame(matrix_id = mid, gse = rec$gse, module = sub("^ME", "", me), module_eigengene = me,
        module_size = sum(colors == sub("^ME", "", me)), soft_power = chosen, n_samples = nrow(datExpr), trait = trait, rho = unname(a$estimate), p_value = a$p.value, stringsAsFactors = FALSE)
    }
  }
}

modules <- rbindlist(module_rows, fill = TRUE); module_traits <- rbindlist(module_trait_rows, fill = TRUE); module_genes <- rbindlist(module_gene_rows, fill = TRUE); correlations <- rbindlist(corr_rows, fill = TRUE)
if (nrow(correlations)) correlations[, p_adj := p.adjust(p_value, method = "BH"), by = matrix_id]
write_tsv(modules, file.path(PROJ, "results", "network", "wgcna_module_manifest.tsv")); write_tsv(module_traits, file.path(PROJ, "results", "network", "wgcna_module_trait_correlations.tsv")); write_tsv(module_genes, file.path(PROJ, "results", "network", "wgcna_module_genes.tsv")); write_tsv(correlations, file.path(PROJ, "results", "network", "candidate_module_correlations.tsv"))

## Summarize candidate-module support without treating a module as a universal
## cross-study object. Require at least 3 independent matrices for consensus.
if (nrow(correlations)) {
  consensus <- correlations[is.finite(rho), .(matrix_n = uniqueN(matrix_id), gse_n = uniqueN(gse), support_n = sum(is.finite(p_adj) & p_adj < cor_fdr),
    median_rho = median(rho, na.rm = TRUE), direction_consistency = max(sum(rho > 0), sum(rho < 0)) / .N), by = .(candidate_id, module_eigengene)]
  consensus[, consensus_supported := matrix_n >= 3 & support_n >= 2 & direction_consistency >= 0.67]
  write_tsv(consensus, file.path(PROJ, "results", "network", "candidate_module_consensus.tsv"))
}

## Rank GPCR candidates by full-transcriptome DE statistics where available.
de <- models[unit == "count" & candidate_id %in% candidate_set & is.finite(stat), .(candidate_id, matrix_id, gse, contrast_type, mutant, time_h, stat, p_adj, log2fc)]
if (nrow(de)) {
  de[, rank_score := sign(stat) * -log10(pmax(p_adj, .Machine$double.xmin))]
  write_tsv(de, file.path(PROJ, "results", "pathway", "candidate_count_de_rankings.tsv"))
}

## Render a compact network plot for strong candidate-module correlations.
if (nrow(correlations)) {
  strong <- correlations[is.finite(p_adj) & p_adj < cor_fdr & abs(rho) >= 0.5]
  if (nrow(strong)) {
    png(file.path(PROJ, "figures", "network", "candidate_module_correlations.png"), width = 1900, height = 1400, res = 180)
    print(ggplot(strong, aes(x = module_eigengene, y = reorder(candidate_id, rho), fill = rho)) + geom_tile() + facet_wrap(~ matrix_id, scales = "free_y") + scale_fill_gradient2() + theme_bw(base_size = 8) + labs(x = "module eigengene", y = "candidate", title = "Strong within-series candidate-module correlations"))
    dev.off()
  }
}

out_files <- c(list.files(file.path(PROJ, "results", "network"), full.names = TRUE), list.files(file.path(PROJ, "results", "pathway"), full.names = TRUE), list.files(file.path(PROJ, "figures", "network"), full.names = TRUE))
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE), file.path(PROJ, "checksums", "step10_network_sha256.tsv"))
cat("eligible_network_matrices:", length(network_mids), "\n"); cat("module_rows:", nrow(modules), "\n"); cat("module_trait_rows:", nrow(module_traits), "\n"); cat("module_gene_rows:", nrow(module_genes), "\n"); cat("candidate_module_correlations:", nrow(correlations), "\n")
if (nrow(correlations)) cat("strong_correlations:", sum(correlations$p_adj < cor_fdr & abs(correlations$rho) >= 0.5, na.rm = TRUE), "\n")
finish_log()
