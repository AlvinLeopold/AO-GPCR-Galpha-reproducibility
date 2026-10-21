#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "05_write_report"
input_paths <- c(
  file.path(PROJ, "metadata", "candidates.tsv"),
  file.path(PROJ, "metadata", "source_file_manifest.tsv"),
  file.path(PROJ, "metadata", "series_summary_atcc.tsv"),
  file.path(PROJ, "results", "expression", "matrix_parse_qc.tsv"),
  file.path(PROJ, "results", "qc", "sample_qc.tsv"),
  file.path(PROJ, "results", "deseq2", "deseq2_status.tsv"),
  file.path(PROJ, "results", "integration", "series_matrix_summary.tsv"),
  file.path(PROJ, "results", "integration", "candidate_cross_study_summary.tsv"),
  file.path(PROJ, "results", "integration", "candidate_priority_shortlist.tsv"),
  file.path(PROJ, "config", "analysis_parameters.tsv")
)
start_log(step, input_paths)

read_tsv <- function(p) if (file.exists(p)) fread(p, sep = "\t", encoding = "UTF-8", fill = TRUE) else data.table()
candidate <- read_tsv(input_paths[1])
manifest <- read_tsv(input_paths[2])
series <- read_tsv(input_paths[3])
parse_qc <- read_tsv(input_paths[4])
sample_qc <- read_tsv(input_paths[5])
deseq_status <- read_tsv(input_paths[6])
matrix_summary <- read_tsv(input_paths[7])
score <- read_tsv(input_paths[8])
short <- read_tsv(input_paths[9])
params <- read_tsv(input_paths[10])

fmt <- function(x, digits = 3) {
  if (!length(x) || is.na(x)) return("NA")
  if (is.numeric(x)) formatC(x, format = "fg", digits = digits, flag = "#") else as.character(x)
}
collapse_or <- function(x, maxn = 8) {
  x <- unique(as.character(x)); x <- x[nzchar(x) & !is.na(x)]
  if (!length(x)) return("未注明")
  if (length(x) > maxn) paste0(paste(x[seq_len(maxn)], collapse = ", "), " 等") else paste(x, collapse = ", ")
}
markdown_escape <- function(x) gsub("\\|", "\\\\|", as.character(x), fixed = FALSE)

report_path <- file.path(PROJ, "results", "integration", "ATCC_24927_GPCR_analysis_report.md")
con <- file(report_path, open = "w", encoding = "UTF-8")
sink(con)
cat("# Arthrobotrys oligospora ATCC 24927 GPCR 表达挖掘报告\n\n")
cat("生成时间（UTC）：", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n\n", sep = "")
cat("## 1. 分析范围与输入\n\n")
cat("本项目严格使用用户在 `Aol_GPCR_determined.xlsx` 中保留的 93 条候选，不重新扩充 ExRe 候选集合。候选包含 ", nrow(candidate), " 条蛋白，其中 ", sum(candidate$tm_count == 7), " 条为 7TM、", sum(candidate$tm_count == 6), " 条为 6TM；3 条 `EVM...` ID 也被保留。\n\n", sep = "")
cat("93 条候选均在 ExRe 的 `Arthrobotrys_oligospora.fasta` 中找到，且序列长度与工作簿完全一致。表达矩阵中的 `AOL_s...`、`gene-AOL_s...` 和省略版本号的基因名按规范化 ID 匹配；原始 ID、匹配方式和序列 SHA256 均保留。\n\n")
cat("数据范围为 12 个明确标注 ATCC 24927 的 GEO 表达系列，共 213 个 GEO 样本。处理后表达文件不等同于原始 FASTQ；本轮优先分析 GEO 已提交的 TPM/FPKM/count 表，未下载约 977 GB 的全部原始测序数据。\n\n")
cat("## 2. GEO 系列与处理设计\n\n")
cat("| GSE | 实际输入矩阵 | 样本数 | 基因型/处理 | 时间 | 备注 |\n|---|---:|---:|---|---|---|\n")
for (g in unique(matrix_summary$gse)) {
  mm <- matrix_summary[gse == g]
  note <- if (g == "GSE253759") "与 GSE215079 高度相似，保留但不视为独立重复" else if (g == "GSE210886") "同一 GSE 中保留 Acy 与 PKA 两个命名块" else if (g == "GSE322740") "Day 5 气生菌丝分化，非时间序列" else ""
  cat("| ", g, " | ", nrow(mm), " | ", sum(mm$n_samples), " | ", markdown_escape(collapse_or(mm$genotypes, 6)),
      " / ", markdown_escape(collapse_or(mm$treatments, 4)), " | ", markdown_escape(collapse_or(mm$times, 6)), " | ", note, " |\n", sep = "")
}
cat("\n详细 GEO 字段见 `metadata/sample_metadata_atcc_series.tsv`；输入文件路径、大小、类型和 SHA256 见 `metadata/source_file_manifest.tsv`。\n\n")
cat("## 3. 解析与质量控制\n\n")
long_path <- file.path(PROJ, "results", "expression", "candidate_expression_long.tsv")
cat("- 解析了 ", nrow(parse_qc), " 个表达/计数矩阵，共 ", sum(parse_qc$expression_columns), " 个矩阵样本列。常规表达表中 3 条 EVM ID 没有可匹配的表达行，按缺失处理而不是填充为 0；GSE322740 另有 2 条 AOL ID 在文件中未出现。\n", sep = "")
cat("- `candidate_expression_long.tsv` 包含 ", format(nrow(read_tsv(long_path)), big.mark = ","), " 行；TPM/FPKM 与 raw counts 分开保留。\n", sep = "")
cat("- 全基因样本 QC 记录 ", nrow(sample_qc), " 个样本-单位组合，结果见 `results/qc/sample_qc.tsv`。\n", sep = "")
cat("- 可用 raw-count 矩阵的 DESeq2 状态：\n\n")
cat("| 矩阵 | 样本数 | 组数 | 保留基因数 | 状态 |\n|---|---:|---:|---:|---|\n")
for (i in seq_len(nrow(deseq_status))) cat("| ", deseq_status$matrix_id[i], " | ", deseq_status$n_samples[i], " | ", deseq_status$n_groups[i], " | ", deseq_status$n_genes_after_filter[i], " | ", deseq_status$status[i], " |\n", sep = "")
cat("\n")
cat("## 4. 统计策略\n\n")
cat("1. TPM/FPKM：在每个矩阵内使用 `log2(value + 0.1)` 计算 WT 诱导、突变体诱导、基因型效应和 genotype × time 交互；重复间使用 Welch t 检验或配对交互检验，BH 校正在矩阵/对比/时间家族内进行。\n")
cat("2. Raw counts：对 GSE301423、GSE301452、GSE303416、GSE306156、GSE322548 的 count 列使用 DESeq2；每个矩阵按 genotype × time 组合建模。\n")
cat("3. 跨研究整合只比较检测广度、效应方向、支持的非潜在重复 GSE 数量和 DESeq2 支持数量，不把 TPM、FPKM、counts 的数值直接合并。GSE253759 在保守计数中排除，但其原始结果完整保留。\n\n")
cat("## 5. 候选优先级\n\n")
cat("优先级分数是可审计的预筛选分数，不是新的显著性检验。Tier 1 要求至少 3 个非潜在重复 GSE 的 WT 诱导支持且方向一致性至少 0.67；Tier 2 包括至少 2 个诱导支持 GSE，或跨 GSE 的基因型/DESeq2 证据；其余为 Tier 3。\n\n")
cat("| Tier | 数量 | 建议用途 |\n|---|---:|---|\n")
tc <- score[, .N, by = priority_tier]
for (i in seq_len(nrow(tc))) cat("| ", tc$priority_tier[i], " | ", tc$N[i], " | ", if (grepl("Tier 1", tc$priority_tier[i])) "优先做独立表达确认、定位和诱导表型" else if (grepl("Tier 2", tc$priority_tier[i])) "按处理条件选择性复核" else "先核对基因模型、检测和重复性", " |\n", sep = "")
cat("\n### Tier 1 候选\n\n")
cat("| 排名 | 候选 | GPCR domain | WT 诱导支持 GSE | 非重复方向一致性 | 基因型支持 GSE | DESeq2 支持 GSE | 方向 | 分数 |\n|---:|---|---|---:|---:|---:|---:|---|---:|\n")
t1 <- score[priority_tier == "Tier 1: multi-GSE supported"]
for (i in seq_len(nrow(t1))) cat("| ", t1$priority_rank[i], " | ", t1$candidate_id[i], " | ", markdown_escape(t1$gpcr_domains[i]), " | ", t1$wt_induction_supported_gse_n_conservative[i], " | ", fmt(t1$wt_gse_direction_consistency_conservative[i]), " | ", t1$genotype_effect_supported_gse_n_conservative[i], " | ", t1$deseq_supported_gse_n_conservative[i], " | ", t1$wt_gse_direction_conservative[i], " | ", fmt(t1$priority_score[i]), " |\n", sep = "")
cat("\n完整 93 条排序见 `results/integration/candidate_priority_shortlist.tsv` 和 `results/integration/candidate_cross_study_summary.tsv`。\n\n")
cat("## 6. 下一步针对性分析计划\n\n")
cat("- **Tier 1（", nrow(t1), "条）**：在 ATCC 24927 同一培养体系中做 qRT-PCR 或 targeted RNA，覆盖 0 h 和线虫/线虫提取物诱导峰值时间；随后做膜定位、过表达/回补和 trap formation 表型。优先从 AOL_s00054g626.1、AOL_s00140g54.6、AOL_s00081g30.1、AOL_s00004g22.1 等开始。\n", sep = "")
cat("- **时间序列**：对 GSE215079/GSE253759（0/6/12/24 h）、GSE245014（0/12/36 h）和 GSE322548（0/12/24 h）绘制逐重复轨迹，区分早期峰、晚期持续响应和缺失响应；先核实 GSE253759 是否为独立或平行提交。\n")
cat("- **通路依赖性**：比较 GSE210886（AoPkaC1/AoAcy）、GSE212735（Aoatg11）、GSE301423/301452（Aoham4/Aoham8/Aochk1）、GSE303416（Aohal4b）、GSE306156（AoOpy2）和 GSE322548（Aohog1）中的 genotype × time interaction。\n")
cat("- **结构联合筛选**：把优先级与 TM 数、GPCR domain、序列长度和 ExRe 哈希联用；6TM 候选保留为独立复核组。\n")
cat("- **原始数据复核**：若需提交级全转录组统计，再针对重点 GSE 下载 FASTQ，锁定参考基因组/注释、文库方向和批次后用 nf-core/rnaseq 或 Salmon 重新定量。\n\n")
cat("## 7. 可复现文件\n\n")
cat("`scripts/01_prepare_inputs.R`、`02_parse_expression.R`、`03_analyze_candidates.R`、`04_integrate_results.R`、`05_write_report.R` 分别负责输入冻结、表达解析、候选分析、跨 GSE 整合和本报告生成；`logs/` 保存命令、起止时间、输入哈希、sessionInfo 和软件版本，`checksums/` 保存校验值。\n\n")
cat("## 8. 重要限制\n\n")
cat("GEO processed expression 表可能来自不同注释版本、定量软件和过滤规则；跨研究分数仅用于候选优先级。表达表缺失标为 `not_in_expression_matrix`，不能解释为生物学零表达。正式功能验证仍需独立培养、独立 RNA 样本和表型/定位实验。\n")
sink(); close(con)

write_tsv(data.frame(path = report_path, sha256 = file_sha256(report_path)), file.path(PROJ, "checksums", "step05_report_sha256.tsv"))
write_tsv(data.frame(artifact = c("report", "candidate_shortlist", "candidate_cross_study_summary"),
                     path = c(report_path, file.path(PROJ, "results", "integration", "candidate_priority_shortlist.tsv"), file.path(PROJ, "results", "integration", "candidate_cross_study_summary.tsv")),
                     sha256 = vapply(c(report_path, file.path(PROJ, "results", "integration", "candidate_priority_shortlist.tsv"), file.path(PROJ, "results", "integration", "candidate_cross_study_summary.tsv")), file_sha256, character(1))),
          file.path(PROJ, "results", "integration", "run_artifacts.tsv"))
cat("report:", report_path, "\n")
cat("tier_counts:\n"); print(score[, .N, by = priority_tier])
finish_log()
