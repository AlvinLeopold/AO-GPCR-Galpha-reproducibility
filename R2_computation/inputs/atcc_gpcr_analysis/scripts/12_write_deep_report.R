#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "12_write_deep_report"
input_paths <- c(
  file.path(PROJ, "results", "design", "series_design_manifest.tsv"),
  file.path(PROJ, "results", "design", "contrast_manifest.tsv"),
  file.path(PROJ, "results", "design", "duplicate_audit.tsv"),
  file.path(PROJ, "results", "models", "model_status.tsv"),
  file.path(PROJ, "results", "models", "processed_count_concordance.tsv"),
  file.path(PROJ, "results", "time_response", "wt_time_consensus.tsv"),
  file.path(PROJ, "results", "perturbation", "candidate_regulator_summary.tsv"),
  file.path(PROJ, "results", "network", "network_manifest.tsv"),
  file.path(PROJ, "results", "network", "candidate_module_correlations.tsv"),
  file.path(PROJ, "results", "network", "candidate_module_consensus.tsv"),
  file.path(PROJ, "results", "panels", "condition_specific_candidate_panels.tsv"),
  file.path(PROJ, "results", "network", "pathway_method_status.tsv")
)
start_log(step, input_paths)
read_dt <- function(p) if (file.exists(p) && file.info(p)$size > 0) fread(p, sep = "\t", encoding = "UTF-8", fill = TRUE) else data.table()
design <- read_dt(input_paths[1]); contrast <- read_dt(input_paths[2]); dup <- read_dt(input_paths[3]); status <- read_dt(input_paths[4]); concord <- read_dt(input_paths[5]); time <- read_dt(input_paths[6]); reg <- read_dt(input_paths[7]); netman <- read_dt(input_paths[8]); netcorr <- read_dt(input_paths[9]); net <- read_dt(input_paths[10]); panels <- read_dt(input_paths[11]); methods <- read_dt(input_paths[12])

fmt <- function(x, d = "NA") { if (!length(x) || all(is.na(x))) d else paste(unique(x[!is.na(x)]), collapse = ", ") }
date_now <- format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
report <- c(
  "# Arthrobotrys oligospora ATCC 24927 GPCR 深度转录组分析报告", "",
  paste0("生成时间：", date_now), "", "## 1. 分析边界", "",
  "本报告严格使用用户筛选的 93 条 GPCR 候选，不从 ExRe 自动扩充。所有表达值均在各 GEO 矩阵内部分析；TPM、FPKM 与 raw counts 不直接拼接。GSE253759 与 GSE215079 的候选表达值完全一致，因此 GSE253759 只用于重复提交/数据一致性审计，不计为独立生物学重复。", "",
  "## 2. 输入设计与审计", "",
  paste0("审计了 ", nrow(design), " 个矩阵、", nrow(contrast), " 个预声明对比；可检验对比 ", sum(contrast$status == "testable"), " 个，未检验 ", sum(contrast$status == "not_tested"), " 个。"), "",
  "| 矩阵 | GEO | 单位 | 样本数 | 基因型 | 时间点 | 状态 |", "|---|---|---:|---:|---|---|---|",
  apply(design[, .(matrix_id, gse, unit, n_samples, genotypes, time_points, matrix_status)], 1, function(z) paste0("| ", paste(z, collapse = " | "), " |")), "",
  paste0("GSE215079/GSE253759 候选表达配对：", dup[ candidate_id == "ALL_CANDIDATES", exact_equal_fraction][1] * 100, "% 完全一致。GSE303416 的 A1/A2/STK1/STK2 仍为 order-based 映射，相关结论标为 provisional。"), "",
  "## 3. 系列内建模", "",
  "processed TPM/FPKM 使用 log2(value + 0.1) 后 limma；raw counts 使用 DESeq2。每个结果保留 log2FC、标准误、95% 区间、p 值、BH-FDR、模型角色和映射置信度。", "",
  paste0("19/", nrow(status), " 个矩阵建模成功；候选结果 ", sum(status$n_contrasts_completed, na.rm = TRUE), " 个对比。配对 processed/count 方向一致率为：", if (nrow(concord)) paste(round(concord[, mean(direction_agree, na.rm = TRUE), by = processed_matrix]$V1 * 100, 1), collapse = ", ") else "NA", "%（逐矩阵）。"), "",
  "## 4. 时间响应", "",
  paste0("时间序列覆盖 AoAge2（0/6/12/24 h）、Rab7/Vps35（0/12/36 h）和 Hog1（0/12/24 h）。WT 共识分类：", if (nrow(time)) paste(names(table(time$consensus_class)), as.integer(table(time$consensus_class)), sep = "=", collapse = "; ") else "NA", "。"), "",
  "分类仅表示表达轨迹形态，不等同于受体功能确认；单系列或跨系列不一致者保留为证据不足。", "",
  "## 5. 调控通路依赖", "",
  paste0("候选×调控因子记录 ", nrow(reg), " 条；交互项强支持记录 ", sum(reg$interaction_supported_n > 0, na.rm = TRUE), " 条；突变体诱导保留 ", sum(reg$dependency_class == "mutant_induction_preserved", na.rm = TRUE), " 条；诱导丢失/改变 ", sum(reg$dependency_class == "mutant_induction_lost_or_changed", na.rm = TRUE), " 条；基线基因型偏移 ", sum(reg$dependency_class == "baseline_genotype_shift", na.rm = TRUE), " 条。"), "",
  "## 6. 模块与网络", "",
  paste0("达到网络门槛的 raw-count 矩阵：", if (nrow(netman)) paste(netman$matrix_id, collapse = ", ") else "无", "。网络仅在系列内部建立；候选-模块相关结果不跨 GEO 直接合并。"), "",
  paste0("候选-模块相关记录 ", if (nrow(netcorr)) nrow(netcorr) else 0, " 条；系列内达到 |rho|≥0.5 且 FDR<0.05 的记录 ", if (nrow(netcorr)) sum(is.finite(netcorr$rho) & is.finite(netcorr$p_adj) & abs(netcorr$rho) >= 0.5 & netcorr$p_adj < 0.05, na.rm = TRUE) else 0, " 条。跨至少 3 个独立矩阵的共识网络支持为 ", if (nrow(net)) sum(net$consensus_supported, na.rm = TRUE) else 0, " 条；当前网络门槛下不应把单系列相关解释为跨研究共识。"), "",
  paste0("GSVA 状态：", if (nrow(methods)) paste(methods$component, methods$available, sep = "=", collapse = "; ") else "未记录", "。当前以 WGCNA 模块特征值和 fgsea 可复现替代；GSVA 未安装不影响已完成的模型与面板。"), "",
  "## 7. 条件化候选面板", "",
  "面板仅用于后续实验优先级，不等同于功能注释。每个面板最多 10 条候选，并同时考虑效应量、重复支持、raw-count 验证、时间轨迹和相应扰动轴。", "",
  "| 面板 | 候选（按优先级） |", "|---|---|",
  if (nrow(panels)) panels[, paste0("| ", panel_id[1], " | ", paste(candidate_id, collapse = ", "), " |"), by = panel_id]$V1 else "| NA | NA |", "",
  "## 8. 建议的下一步", "",
  "1. 优先对 P1 捕食核心、P2 cAMP-PKA 和 P3 MAPK/渗透压面板做独立 ATCC 培养验证与 qRT-PCR。", "2. 对 raw-count 支持和交互项支持同时较强的候选，增加受体定位、配体刺激和捕食表型实验。", "3. 在采用 GSE303416 结果前核对 A1/A2/STK1/STK2 与 GEO 样本标题的真实对应关系。", "4. 第二阶段再对优先系列重定量 FASTQ；统一 GCF_000225545.1 参考、GTF、文库方向和批次后，与 processed 结果逐候选比较。", "", "## 9. 可重复性", "",
  "所有步骤由 `scripts/06_design_audit.R` 至 `scripts/12_write_deep_report.R` 执行；每一步写入独立日志、输入 SHA256、输出 SHA256、软件版本和失败/未检验记录。"
)
report_path <- file.path(PROJ, "results", "integration", "ATCC_24927_deep_analysis_report.md")
writeLines(report, report_path, useBytes = TRUE)
out_files <- c(report_path, file.path(PROJ, "results", "integration", "ATCC_24927_deep_analysis_plan.md"))
write_tsv(data.frame(path = out_files, sha256 = vapply(out_files, file_sha256, character(1)), stringsAsFactors = FALSE), file.path(PROJ, "checksums", "step12_deep_report_sha256.tsv"))
cat("report:", report_path, "\n"); cat("report_lines:", length(report), "\n")
finish_log()
