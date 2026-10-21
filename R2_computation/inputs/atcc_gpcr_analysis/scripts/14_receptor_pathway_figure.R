#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(readxl)
  library(scales)
})

source(file.path(getwd(), "atcc_gpcr_analysis", "scripts", "common.R"), local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript <script> <orthogroup_mapping.tsv> <prior_study_supplement.xlsx>")
}
orthogroup_path <- normalizePath(args[1], mustWork = TRUE)
prior_supplement_path <- normalizePath(args[2], mustWork = TRUE)

model_path <- file.path(PROJ, "results", "models", "candidate_series_model_results.tsv")
candidate_path <- file.path(PROJ, "metadata", "candidates.tsv")
source_dir <- file.path(PROJ, "results", "figure_source_data")
figure_dir <- file.path(PROJ, "figures", "manuscript")
safe_mkdir(source_dir)
safe_mkdir(figure_dir)

targets <- c(
  "AOL_s00007g309.1", "AOL_s00054g150.1", "AOL_s00176g42.1",
  "AOL_s00054g200.1", "AOL_s00097g69.1", "AOL_s00110g16.1",
  "AOL_s00140g54.6", "AOL_s00173g37.1", "AOL_s00173g137.1"
)
focus_targets <- targets[1:3]

backgrounds <- data.table(
  matrix_id = c(
    "GSE210886_Acy", "GSE210886_PKA", "GSE306156_TPM", "GSE322548_TPM",
    "GSE245014_Rab7", "GSE245014_Vps35", "GSE212735", "GSE213542",
    "GSE215079", "GSE301423_TPM", "GSE301423_TPM", "GSE301452_FPKM",
    "GSE303416_FPKM"
  ),
  mutant = c(
    "delta_AoAcy", "delta_AoPkaC1", "delta_AoOpy2", "delta_Aohog1",
    "delta_AoRab7", "delta_AoVps35", "delta_Aoatg11", "delta_AobZip129",
    "delta_Aoage2", "delta_Aoham4", "delta_Aoham8", "delta_Aochk1",
    "delta_Aohal4b"
  ),
  regulator = c(
    "AoAcy", "AoPkaC1", "AoOpy2", "AoHog1", "AoRab7", "AoVps35",
    "AoAtg11", "AobZip129", "AoAge2", "AoHam4", "AoHam8", "AoChk1",
    "AoHal4b"
  ),
  pathway_group = c(
    "cAMP-PKA", "cAMP-PKA", "Osmotic MAPK", "Osmotic MAPK",
    "Vacuolar trafficking", "Vacuolar trafficking", rep("Other regulatory backgrounds", 7)
  ),
  background_order = seq_len(13)
)

models <- fread(model_path, sep = "\t", encoding = "UTF-8")
candidates <- fread(candidate_path, sep = "\t", encoding = "UTF-8")
selected <- merge(
  models[
    candidate_id %chin% targets & contrast_type == "genotype_effect" & unit != "count",
    .(candidate_id, matrix_id, gse, unit, method, contrast_type, mutant, time_h,
      log2fc, se, ci_low, ci_high, p_value, p_adj, mapping_confidence, model_role)
  ],
  backgrounds,
  by = c("matrix_id", "mutant"),
  all = FALSE
)

expected_pairs <- CJ(candidate_id = targets, background_order = seq_len(nrow(backgrounds)))
observed_pairs <- unique(selected[, .(candidate_id, background_order)])
if (nrow(observed_pairs) != nrow(expected_pairs) ||
    nrow(merge(expected_pairs, observed_pairs, by = c("candidate_id", "background_order"))) != nrow(expected_pairs)) {
  stop("The 9 x 13 receptor-background grid is incomplete.")
}

context <- selected[, .(
  matrix_id = first(matrix_id),
  gse = first(gse),
  unit = first(unit),
  method = first(method),
  pathway_group = first(pathway_group),
  regulator = first(regulator),
  timepoint_n = .N,
  timepoints_h = paste(sort(unique(time_h)), collapse = ";"),
  median_log2fc = median(log2fc, na.rm = TRUE),
  mean_log2fc = mean(log2fc, na.rm = TRUE),
  min_log2fc = min(log2fc, na.rm = TRUE),
  max_log2fc = max(log2fc, na.rm = TRUE),
  max_abs_log2fc = max(abs(log2fc), na.rm = TRUE),
  supported_time_n = sum(is.finite(p_adj) & p_adj < 0.05),
  any_fdr = any(is.finite(p_adj) & p_adj < 0.05),
  strong_time_n = sum(is.finite(p_adj) & p_adj < 0.05 & abs(log2fc) >= 1),
  any_strong = any(is.finite(p_adj) & p_adj < 0.05 & abs(log2fc) >= 1)
), by = .(candidate_id, background_order)]

context[, `:=`(
  background_mean_abs_median = mean(abs(median_log2fc)),
  background_fdr_receptor_n = sum(any_fdr),
  background_strong_receptor_n = sum(any_strong),
  background_up_n = sum(median_log2fc > 0),
  background_down_n = sum(median_log2fc < 0)
), by = background_order]

structural <- candidates[
  sequence_id %chin% targets,
  .(candidate_id = sequence_id, tm_count, length_workbook, gpcr_domain_count, gpcr_domains)
]
context <- merge(context, structural, by = "candidate_id", all.x = TRUE)
context[, display_id := sub("\\.[0-9]+$", "", candidate_id)]
context[, focus_gene := candidate_id %chin% focus_targets]
context[, focus_pair :=
  (candidate_id == "AOL_s00007g309.1" & regulator %chin% c("AoHog1", "AoRab7", "AoVps35")) |
  (candidate_id == "AOL_s00054g150.1" & regulator == "AoHog1") |
  (candidate_id == "AOL_s00176g42.1" & regulator %chin% c("AoAcy", "AoVps35"))
]
context[, effect_definition := paste(
  "Median across matched-time mutant-versus-WT genotype effects;",
  "background ranking is mean absolute receptor-level median across the nine receptors"
)]
context[, target_order := match(candidate_id, targets)]
setorder(context, background_order, target_order)

write_tsv(
  context,
  file.path(source_dir, "Fig4_receptor_pathway_context.tsv")
)

key_processed <- copy(selected)
key_processed[, evidence_role := "primary_processed_expression_for_13_background_summary"]
key_processed[, display_id := sub("\\.[0-9]+$", "", candidate_id)]
key_processed[, focus_pair :=
  (candidate_id == "AOL_s00007g309.1" & regulator %chin% c("AoHog1", "AoRab7", "AoVps35")) |
  (candidate_id == "AOL_s00054g150.1" & regulator == "AoHog1") |
  (candidate_id == "AOL_s00176g42.1" & regulator %chin% c("AoAcy", "AoVps35"))
]

hog_counts <- models[
  candidate_id %chin% targets & matrix_id == "GSE322548_counts" &
    contrast_type == "genotype_effect" & unit == "count",
  .(candidate_id, matrix_id, gse, unit, method, contrast_type, mutant, time_h,
    log2fc, se, ci_low, ci_high, p_value, p_adj, mapping_confidence, model_role)
]
hog_counts[, `:=`(
  regulator = "AoHog1",
  pathway_group = "Osmotic MAPK",
  background_order = 4L,
  evidence_role = "confirmatory_raw_count_sensitivity",
  display_id = sub("\\.[0-9]+$", "", candidate_id),
  focus_pair = candidate_id %chin% c("AOL_s00007g309.1", "AOL_s00054g150.1")
)]

key_columns <- c(
  "candidate_id", "display_id", "regulator", "pathway_group", "background_order",
  "matrix_id", "gse", "unit", "method", "contrast_type", "mutant", "time_h",
  "log2fc", "se", "ci_low", "ci_high", "p_value", "p_adj", "mapping_confidence",
  "model_role", "evidence_role", "focus_pair"
)
key_evidence <- rbindlist(list(key_processed[, ..key_columns], hog_counts[, ..key_columns]), fill = TRUE)
setorder(key_evidence, evidence_role, background_order, candidate_id, time_h)
write_tsv(
  key_evidence,
  file.path(source_dir, "Fig4_receptor_key_evidence.tsv")
)

orthogroups <- fread(orthogroup_path, sep = "\t", encoding = "UTF-8")
bridge <- orthogroups[
  vapply(AOL_protein_ids, function(z) any(targets %chin% trimws(strsplit(z, ",", fixed = TRUE)[[1]])), logical(1))
]
bridge_long <- rbindlist(lapply(seq_len(nrow(bridge)), function(i) {
  aol_ids <- trimws(strsplit(bridge$AOL_protein_ids[i], ",", fixed = TRUE)[[1]])
  aol_id <- intersect(aol_ids, targets)[1]
  twf_ids <- trimws(strsplit(bridge$TWF_protein_ids[i], ",", fixed = TRUE)[[1]])
  data.table(
    candidate_id = aol_id,
    orthogroup = bridge$orthogroup[i],
    relation_type = bridge$relation_type[i],
    explicit_orthologues_call = bridge$explicit_orthologues_call[i],
    twf_id = twf_ids
  )
}))

prior_annotation_raw <- as.data.table(read_excel(
  prior_supplement_path, sheet = "Suppl.Table 1", skip = 1, .name_repair = "minimal"
))
prior_annotation <- prior_annotation_raw[
  grepl("^EYR41_", as.character(prior_annotation_raw[[1]])),
  .(twf_id = as.character(prior_annotation_raw[[1]]),
    prior_classification = as.character(prior_annotation_raw[[3]]))
]
prior_response_raw <- as.data.table(read_excel(
  prior_supplement_path, sheet = "Suppl.Table 2", skip = 1, .name_repair = "minimal"
))
prior_response <- prior_response_raw[
  grepl("^EYR41_", as.character(prior_response_raw[[1]])),
  .(twf_id = as.character(prior_response_raw[[1]]),
    prior_tpm_0h = as.numeric(prior_response_raw[[2]]),
    prior_tpm_10h = as.numeric(prior_response_raw[[3]]))
]
# Ratios are defined only for strictly positive mean-TPM values.
if (any(prior_response$prior_tpm_0h <= 0 | prior_response$prior_tpm_10h <= 0, na.rm = TRUE)) {
  stop("Prior-study mean TPM values must be positive before calculating log2 ratios.")
}
prior_response[, prior_log2fc_10h_vs_0h := log2(prior_tpm_10h / prior_tpm_0h)]

prior_bridge <- merge(bridge_long, prior_annotation, by = "twf_id", all.x = TRUE)
prior_bridge <- merge(prior_bridge, prior_response, by = "twf_id", all.x = TRUE)
prior_bridge[, display_id := sub("\\.[0-9]+$", "", candidate_id)]
prior_bridge[, mapping_evidence := fifelse(
  explicit_orthologues_call == "yes", "Explicit orthologue call", "Orthogroup-only bridge"
)]
prior_bridge[, response_table_status := fifelse(
  is.finite(prior_log2fc_10h_vs_0h), "Mean-TPM estimate available", "Not present in the supplied response table"
)]
prior_bridge[, target_order := match(candidate_id, targets)]
setorder(prior_bridge, target_order, twf_id)
write_tsv(
  prior_bridge,
  file.path(source_dir, "Fig4_prior_study_bridge.tsv")
)

font_family <- "Helvetica"
palette <- list(
  ink = "#202124",
  muted = "#667078",
  light = "#EEF1F2",
  up = "#C95D3A",
  down = "#356A92",
  mapk = "#D87C4A",
  camp = "#4379A8",
  vac = "#4F8A68",
  accent = "#D6A11D",
  association = "#0E6F70"
)

theme_set(theme_classic(base_size = 6.4, base_family = font_family))

box_layer <- function(x, y, width, height, fill, colour = palette$ink) {
  annotate(
    "rect", xmin = x - width / 2, xmax = x + width / 2,
    ymin = y - height / 2, ymax = y + height / 2,
    fill = fill, colour = colour, linewidth = 0.35
  )
}
arrow_layer <- function(x, xend, y, yend, colour = palette$ink, linetype = "solid") {
  annotate(
    "segment", x = x, xend = xend, y = y, yend = yend,
    colour = colour, linewidth = 0.5, linetype = linetype,
    arrow = arrow(length = unit(1.5, "mm"), type = "closed")
  )
}

p_mechanism <- ggplot() +
  coord_cartesian(xlim = c(0, 13.4), ylim = c(0, 8), clip = "off") +
  box_layer(1.15, 6.0, 1.75, 1.25, "#F2F3F4") +
  box_layer(4.0, 6.0, 2.45, 1.6, "#F2F3F4") +
  box_layer(6.75, 6.0, 1.65, 1.25, "#F2F3F4") +
  box_layer(9.15, 6.0, 1.75, 1.25, "#EAF1F7", palette$camp) +
  box_layer(12.0, 6.0, 1.95, 1.25, "#F2F3F4") +
  arrow_layer(2.03, 2.92, 6.48, 6.48) +
  arrow_layer(5.23, 5.92, 6.48, 6.48) +
  arrow_layer(7.58, 8.28, 6.48, 6.48) +
  arrow_layer(10.03, 11.03, 6.48, 6.48) +
  box_layer(5.0, 4.00, 3.45, 1.35, "#FFF4DE", palette$accent) +
  box_layer(1.70, 1.65, 3.25, 1.95, "#FAEEE7", palette$mapk) +
  box_layer(9.65, 1.65, 3.35, 1.95, "#EAF3ED", palette$vac) +
  annotate(
    "segment", x = 3.10, xend = 3.65, y = 2.48, yend = 3.42,
    colour = palette$association, linewidth = 0.55, linetype = "22"
  ) +
  annotate(
    "segment", x = 3.42, xend = 3.92, y = 3.38, yend = 3.38,
    colour = palette$association, linewidth = 0.7
  ) +
  annotate(
    "segment", x = 8.25, xend = 6.72, y = 2.48, yend = 3.48,
    colour = palette$association, linewidth = 0.55, linetype = "22"
  ) +
  annotate(
    "segment", x = 8.45, xend = 6.78, y = 5.38, yend = 4.62,
    colour = palette$association, linewidth = 0.55, linetype = "22"
  ) +
  annotate("text", x = 1.15, y = 6.0, label = "Nematode\ncues", size = 1.92, family = font_family) +
  annotate("text", x = 4.0, y = 6.18, label = "Reported GPCRs", size = 2.0, fontface = "bold", family = font_family) +
  annotate("text", x = 4.0, y = 5.70, label = "Gpr2, Gpr3 and Gin3", size = 1.72, family = font_family, colour = palette$muted) +
  annotate("text", x = 6.75, y = 6.0, label = "G protein\nGpa2", size = 1.92, family = font_family) +
  annotate("text", x = 9.15, y = 6.0, label = "cAMP-PKA", size = 2.05, fontface = "bold", family = font_family, colour = palette$camp) +
  annotate("text", x = 12.0, y = 6.0, label = "Trap development\nPredatory transition", size = 1.72, family = font_family) +
  annotate("text", x = 5.0, y = 4.22, label = "9 ATCC 24927 GPCR candidates", size = 1.92, fontface = "bold", family = font_family) +
  annotate("text", x = 5.0, y = 3.80, label = "7-TM; Pth11-like domains", size = 1.70, family = font_family, colour = palette$muted) +
  annotate("text", x = 4.95, y = 4.88, label = "inferred negative\ntranscriptional control", size = 1.68, family = font_family, colour = palette$association, hjust = 1) +
  annotate("text", x = 1.70, y = 2.22, label = "AoOpy2 / AoHog1", size = 1.96, fontface = "bold", family = font_family, colour = palette$mapk) +
  annotate("text", x = 1.70, y = 1.78, label = "Osmotic MAPK context", size = 1.72, family = font_family, colour = palette$muted) +
  annotate("text", x = 1.70, y = 1.30, label = "9/9 FDR; 6 up / 3 down\nmean |median| = 1.95", size = 1.52, family = font_family, colour = palette$mapk) +
  annotate("text", x = 9.65, y = 2.22, label = "AoRab7 / AoVps35", size = 1.96, fontface = "bold", family = font_family, colour = palette$vac) +
  annotate("text", x = 9.65, y = 1.78, label = "Vacuolar trafficking context", size = 1.72, family = font_family, colour = palette$muted) +
  annotate("text", x = 9.65, y = 1.30, label = "Receptor-specific bidirectional effects", size = 1.52, family = font_family, colour = palette$vac) +
  annotate("text", x = 5.65, y = 0.45, label = "Solid: previously reported route", size = 1.80, family = font_family, colour = palette$muted, hjust = 1) +
  annotate("segment", x = 5.8, xend = 6.25, y = 0.45, yend = 0.45, colour = palette$association, linewidth = 0.55, linetype = "22") +
  annotate("text", x = 6.4, y = 0.45, label = "Dashed: transcriptome association / working model", size = 1.80, family = font_family, colour = palette$association, hjust = 0) +
  labs(title = "Reported prey-sensing route and receptor-centred working model") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(3, 2, 1, 2)
  )

prior_plot <- prior_bridge[is.finite(prior_log2fc_10h_vs_0h)]
prior_plot <- prior_plot[order(prior_log2fc_10h_vs_0h)]
prior_plot[, display_id := factor(display_id, levels = unique(display_id))]
p_prior <- ggplot(prior_plot, aes(x = prior_log2fc_10h_vs_0h, y = display_id)) +
  geom_vline(xintercept = 0, colour = "#9AA0A6", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = prior_log2fc_10h_vs_0h, yend = display_id), colour = "#C7CDD1", linewidth = 0.55) +
  geom_point(aes(fill = prior_log2fc_10h_vs_0h, shape = mapping_evidence), colour = palette$ink, size = 2.2, stroke = 0.55) +
  geom_text(aes(label = sprintf("%+.1f", prior_log2fc_10h_vs_0h)), hjust = -0.35, size = 2.0, family = font_family) +
  scale_fill_gradient2(low = palette$down, mid = "white", high = palette$up, midpoint = 0, limits = c(-6.5, 6.5), oob = squish, guide = "none") +
  scale_shape_manual(values = c("Explicit orthologue call" = 21, "Orthogroup-only bridge" = 24)) +
  scale_x_continuous(limits = c(-1.0, 7.25), breaks = c(0, 2, 4, 6), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "TWF154 exposure-linked counterparts",
    subtitle = "Reported mean TPM, 10 h versus 0 h",
    x = "log2 fold change", y = NULL, shape = NULL
  ) +
  theme_classic(base_size = 6.2, base_family = font_family) +
  theme(
    plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 5.8, colour = palette$muted, margin = margin(b = 2)),
    axis.text.y = element_text(size = 5.5),
    axis.text.x = element_text(size = 5.5),
    axis.title.x = element_text(size = 5.8),
    legend.position = "none",
    plot.margin = margin(3, 2, 1, 2)
  )

heat <- copy(context)
heat[, display_id := factor(display_id, levels = rev(sub("\\.[0-9]+$", "", targets)))]
heat[, regulator := factor(regulator, levels = backgrounds$regulator)]
heat[, pathway_group := factor(
  pathway_group,
  levels = c("cAMP-PKA", "Osmotic MAPK", "Vacuolar trafficking", "Other regulatory backgrounds")
)]
levels(heat$pathway_group) <- c("cAMP-PKA", "HOG-MAPK", "Vacuolar", "Other regulatory backgrounds")
heat[, tile_label := sprintf("%+.1f%s", median_log2fc, ifelse(any_fdr, "*", ""))]
focus_heat <- heat[focus_pair == TRUE]

p_heat <- ggplot(heat, aes(x = regulator, y = display_id, fill = median_log2fc)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_tile(data = focus_heat, fill = NA, colour = palette$accent, linewidth = 0.8) +
  geom_text(aes(label = tile_label), size = 1.80, family = font_family, colour = "#18191A") +
  facet_grid(cols = vars(pathway_group), scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = palette$down, mid = "#FAFAFA", high = palette$up,
    midpoint = 0, limits = c(-5, 5), oob = squish,
    breaks = c(-4, -2, 0, 2, 4), name = "Median\nlog2FC"
  ) +
  labs(
    title = "Nine-receptor response landscape across 13 genetic backgrounds",
    subtitle = "Matched-time mutant versus WT; cell = cross-time median; * at least one FDR < 0.05",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 6.2, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.1, "mm"),
    strip.background = element_rect(fill = "#EEF1F2", colour = NA),
    strip.text.x = element_text(size = 5.6, face = "bold", margin = margin(1.5, 1, 1.5, 1)),
    axis.text.x = element_text(size = 5.3, angle = 48, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5.6, colour = palette$ink),
    legend.position = "bottom",
    legend.key.width = unit(17, "mm"),
    legend.key.height = unit(2.1, "mm"),
    legend.title = element_text(size = 5.6),
    legend.text = element_text(size = 5.3),
    plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 5.7, colour = palette$muted, margin = margin(b = 3)),
    plot.margin = margin(3, 2, 2, 2)
  )

focus_order <- c(
  "AOL_s00007g309 | AoRab7",
  "AOL_s00007g309 | AoVps35",
  "AOL_s00176g42 | AoVps35",
  "AOL_s00176g42 | AoAcy",
  "AOL_s00054g150 | AoHog1",
  "AOL_s00007g309 | AoHog1"
)
focus <- context[focus_pair == TRUE]
focus[, pair_label := paste(display_id, regulator, sep = " | ")]
focus[, pair_label := factor(pair_label, levels = focus_order)]
focus[, direction := fifelse(median_log2fc >= 0, "Up", "Down")]
focus[, label_x := fifelse(
  median_log2fc >= 0,
  pmin(5.95, max_log2fc + 0.25),
  pmax(-4.35, min_log2fc - 0.25)
)]
focus[, label_hjust := fifelse(
  median_log2fc >= 0 & label_x >= 5.9, 1,
  fifelse(median_log2fc < 0 & label_x <= -4.3, 0, fifelse(median_log2fc >= 0, 0, 1))
)]

p_focus <- ggplot(focus, aes(y = pair_label, x = median_log2fc)) +
  geom_vline(xintercept = 0, colour = "#9AA0A6", linewidth = 0.35) +
  geom_segment(aes(x = min_log2fc, xend = max_log2fc, yend = pair_label), colour = "#737B80", linewidth = 0.8, lineend = "round") +
  geom_point(aes(fill = direction), shape = 21, colour = palette$ink, size = 2.35, stroke = 0.5) +
  geom_text(aes(x = label_x, label = sprintf("%+.1f", median_log2fc), hjust = label_hjust), size = 2.0, family = font_family) +
  scale_fill_manual(values = c("Up" = palette$up, "Down" = palette$down), guide = "none") +
  scale_x_continuous(limits = c(-4.6, 6.2), breaks = c(-4, -2, 0, 2, 4, 6), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Focal pathway associations",
    subtitle = "Point: median; line: range across sampled times",
    x = "Genotype effect (log2FC)", y = NULL
  ) +
  theme_classic(base_size = 6.2, base_family = font_family) +
  theme(
    plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 5.7, colour = palette$muted, margin = margin(b = 2)),
    axis.text.y = element_text(size = 5.25),
    axis.text.x = element_text(size = 5.4),
    axis.title.x = element_text(size = 5.8),
    plot.margin = margin(3, 4, 2, 2)
  )

top_row <- p_mechanism + p_prior + plot_layout(widths = c(1.75, 1))
bottom_row <- p_heat + p_focus + plot_layout(widths = c(1.78, 1))
figure <- top_row / bottom_row +
  plot_layout(heights = c(0.82, 1.45)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 8, face = "bold", family = font_family),
      plot.tag.position = c(0, 1),
      plot.margin = margin(2, 2, 2, 2)
    )
  )

base_out <- file.path(figure_dir, "Fig4_receptor_pathway_context")
width_mm = 183
height_mm = 150
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(base_out, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()

grDevices::cairo_pdf(paste0(base_out, ".pdf"), width = width_in, height = height_in, family = "Helvetica")
print(figure)
dev.off()

ragg::agg_tiff(
  paste0(base_out, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure)
dev.off()

ragg::agg_png(
  paste0(base_out, ".png"), width = width_in, height = height_in,
  units = "in", res = 300
)
print(figure)
dev.off()

cat("Figure:", paste0(base_out, ".{svg,pdf,tiff,png}"), "\n")
cat("Context rows:", nrow(context), "\n")
cat("Key evidence rows:", nrow(key_evidence), "\n")
cat("Prior bridge rows:", nrow(prior_bridge), "\n")
cat("AoHog1 mean absolute receptor median:", sprintf("%.3f", context[regulator == "AoHog1", first(background_mean_abs_median)]), "\n")
cat("AoAcy mean absolute receptor median:", sprintf("%.3f", context[regulator == "AoAcy", first(background_mean_abs_median)]), "\n")
