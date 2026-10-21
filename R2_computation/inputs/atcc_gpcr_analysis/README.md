# ATCC 24927 GPCR transcriptome analysis

本项目针对 *Arthrobotrys oligospora* ATCC 24927，使用同级目录 `../Aol_GPCR_determined.xlsx` 中筛选后的 93 条候选 GPCR。候选序列从 `../reference/Arthrobotrys_oligospora.fasta` 提取；表达数据来自 `../geo_ao_transcriptomes/catalog/raw/supplementary/` 下 12 个 GEO 系列的已下载 processed expression/count 文件。这些是 GEO 提供的表达/计数矩阵，**不是原始 FASTQ reads**。

## Scope

- 93 条候选：90 条 7TM、3 条 6TM；3 条 `EVM...` ID 保留。
- 12 个 ATCC GEO 系列，213 个 GEO 样本。
- 15 个来源文件，对应 19 个实际表达/计数矩阵。
- TPM、FPKM 和 raw counts 分开分析，不跨单位直接合并数值。
- GSE253759 与 GSE215079 的设计/列名高度相似；保留全部数据，但保守跨 GSE 支持计数排除 GSE253759。

## Reproduce

在本项目目录执行（须预先安装 R 4.4.1 及所需 CRAN/Bioconductor 包；原工作站的 `sessionInfo()` 见 `logs/`）：

```bash
./run_all.sh
```

脚本按顺序运行 `scripts/01_prepare_inputs.R` 到 `scripts/13_nature_figures.R`。其中 `06`--`12` 分别执行设计审计、系列内 limma/DESeq2、时间响应、扰动共识、WGCNA/网络、条件化面板和深度报告；`13` 使用 R（ggplot2/patchwork/ragg/svglite）生成 Nature 风格主图、可编辑 SVG/PDF、600 dpi TIFF 和预览 PNG。每一步生成 `logs/<step>.log`，记录命令、UTC 起止时间、输入文件 SHA256、R `sessionInfo()` 和软件版本。`checksums/` 保存输入和主要输出校验值。

## Main outputs

- `metadata/candidates.tsv`：93 条候选、工作簿字段、ExRe 序列长度和 SHA256。
- `metadata/candidate_sequences.fasta`：严格的 93 条候选 FASTA 子集。
- `metadata/sample_metadata_atcc_series.tsv`：12 个 GSE 的完整 GEO 样本设计。
- `metadata/source_file_manifest.tsv`：每个输入文件的类型、大小和 SHA256。
- `results/expression/candidate_expression_long.tsv`：候选-样本长表，含原始 ID、匹配方式和单位。
- `results/qc/`：全基因样本 QC 和候选重复离散程度。
- `results/effects/candidate_effects.tsv`：WT 诱导、突变体诱导、基因型效应和交互效应。
- `results/deseq2/deseq2_candidate_results.tsv`：count-bearing 矩阵的 DESeq2 候选结果。
- `results/integration/candidate_cross_study_summary.tsv`：93 条候选跨 GSE 证据总表。
- `results/integration/candidate_priority_shortlist.tsv`：预筛选分层和排序。
- `results/integration/ATCC_24927_GPCR_analysis_report.md`：中文总结报告和后续实验计划。
- `results/integration/ATCC_24927_deep_analysis_report.md`：深度分析设计审计、系列模型、时间响应、扰动共识、网络和条件化面板总结。
- `results/design/`：系列设计、对比清单、样本映射、重复提交和 GSM 重叠审计。
- `results/models/`：系列内 limma/DESeq2 结果、PCA、processed/count 一致性和 raw-count 全基因结果。
- `results/time_response/`：时间轨迹、峰值/AUC、WT/突变体/交互分型和共识分类。
- `results/perturbation/`：调控因子证据、配对扰动比较和候选通路指纹。
- `results/network/`：系列内 WGCNA 模块、候选-模块相关和方法状态；GSVA 缺失已显式记录。
- `results/panels/`：六类条件化候选面板、分数和面板成员关系。
- `figures/`：样本信号、基线表达、WT 诱导方向和优先级图；`figures/nature/` 为三组提交级组合图。
- `results/figure_source_data/`：每个图 panel 的逐行源数据；`results/figure_legends.md` 和 `results/figure_contract.md` 记录图注、数据范围和解释边界；`results/figure_QA.md` 记录尺寸、分辨率和视觉检查。
- `results/figure_descriptions.md`：三组主图的英文 title、figure description/legend 和 Results 段落，可直接作为论文初稿使用。

## Interpretation

候选优先级是可重复的预筛选排序，不等同于功能确认。processed expression 文件可能来自不同注释版本、定量软件和过滤规则；表达表中缺失的候选标为 `not_in_expression_matrix`，不能解释为生物学零表达。提交级全转录组结论仍应对重点 GSE 下载 FASTQ，在锁定参考基因组、注释、文库方向和批次后重新定量。
