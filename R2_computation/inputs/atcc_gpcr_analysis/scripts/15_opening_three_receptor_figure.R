#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(grid)
  library(scales)
})

source(file.path(getwd(), "atcc_gpcr_analysis", "scripts", "common.R"), local = TRUE)

figure_dir <- file.path(PROJ, "figures", "manuscript")
source_dir <- file.path(PROJ, "results", "figure_source_data")
safe_mkdir(figure_dir)
safe_mkdir(source_dir)

source_path <- file.path(source_dir, "Fig_opening_three_receptors_context.tsv")
selection_path <- file.path(source_dir, "Fig_opening_three_receptor_selection.tsv")
selection_criterion_path <- file.path(source_dir, "Fig_opening_three_receptor_selection_candidates.tsv")

raw_context <- fread(source_path, encoding = "UTF-8")
raw_selection <- fread(selection_criterion_path, encoding = "UTF-8")

focal_ids <- c("AOL_s00007g309", "AOL_s00054g150", "AOL_s00176g42")
focal_transcripts <- paste0(focal_ids, ".1")
context <- raw_context[
  receptor %chin% focal_ids &
    regulator_background %chin% c("ΔAoHog1", "ΔAoOpy2", "ΔAoAcy", "ΔAoPkaC1", "ΔAoRab7", "ΔAoVps35")
]
context[, display_id := factor(receptor, levels = rev(focal_ids))]
context[, regulator_label := factor(
  regulator_background,
  levels = c("ΔAoHog1", "ΔAoOpy2", "ΔAoAcy", "ΔAoPkaC1", "ΔAoRab7", "ΔAoVps35"),
  labels = c("AoHog1", "AoOpy2", "AoAcy", "AoPkaC1", "AoRab7", "AoVps35")
)]
context[, module := factor(
  pathway_module,
  levels = c("osmosensing MAPK", "cAMP–PKA", "vacuolar trafficking"),
  labels = c("Osmosensing MAPK", "cAMP-PKA", "Vacuolar trafficking")
)]
context[, supported := is.finite(n_supported_BHq005) & n_supported_BHq005 > 0]
context[, tile_label := sprintf("%+.1f%s", median_log2FC_mutant_vs_WT, ifelse(supported, "*", ""))]
context[, regulator_label := factor(regulator_label, levels = c("AoHog1", "AoOpy2", "AoAcy", "AoPkaC1", "AoRab7", "AoVps35"))]

if (nrow(context) != 18L) stop("Expected 3 x 6 focal receptor-background cells.")

selection <- raw_selection
selection_counts <- data.table(
  stage = factor(c("Evaluable GPCRs", "Retained bridge to prior 83", "Prior 10 h induced", "Both MAPK and cAMP-PKA", "Uncharacterized after prior work", "Focal receptors"),
                 levels = c("Evaluable GPCRs", "Retained bridge to prior 83", "Prior 10 h induced", "Both MAPK and cAMP-PKA", "Uncharacterized after prior work", "Focal receptors")),
  n = c(90L, 72L, 29L, 16L, 13L, 3L),
  detail = c("all candidates with modelled expression", "eligibility bridge used for the fixed filter", "within the prior 32-gene response set", "at least one BH-supported contrast in each module", "after excluding prior named/deleted GIN2, GIN6 and GIN9", "largest-effect members within the final 13")
)

exposure <- data.table(
  receptor = focal_ids,
  display_id = factor(focal_ids, levels = rev(focal_ids)),
  twf_id = c("EYR41_011201", "EYR41_002638", "EYR41_010125"),
  prior_tpm_0h = c(13.3650919, 0.02126184, 16.9519923),
  prior_tpm_10h = c(25.953196, 1.55281683, 13.534976),
  cross_geo_support = c(2L, 6L, 4L),
  cross_geo_total = 11L
)
exposure[, prior_log2ratio := log2(prior_tpm_10h / prior_tpm_0h)]
exposure[, direction := ifelse(prior_log2ratio >= 0, "Up", "Down")]

raw_context[, source_table := "provided focal receptor source table"]
exposure_source <- copy(exposure)
exposure_source[, source_table := "Kuo et al. 2024 mean TPM and ATCC cross-GEO support"]
fwrite(context, file.path(source_dir, "Fig_opening_three_receptors_context.tsv"), sep = "\t", na = "NA")
fwrite(selection_counts, file.path(source_dir, "Fig_opening_three_receptor_selection.tsv"), sep = "\t", na = "NA")
fwrite(raw_selection, file.path(source_dir, "Fig_opening_three_receptor_selection_candidates.tsv"), sep = "\t", na = "NA")
fwrite(exposure_source, file.path(source_dir, "Fig_opening_three_receptor_exposure.tsv"), sep = "\t", na = "NA")

font_family <- "Helvetica"
pal <- list(
  ink = "#202124", muted = "#68737A", light = "#F1F3F4", pale_blue = "#E9F0F6",
  up = "#C43D4B", down = "#356A92", mapk = "#B65A3B", camp = "#3F78A9",
  vac = "#4B856B", teal = "#167C7B", gold = "#C7922D", grey = "#B8C0C5"
)
theme_set(theme_classic(base_size = 6.6, base_family = font_family))

arrow_segment <- function(x, xend, y, yend, colour = pal$ink, linetype = "solid", linewidth = 0.55) {
  annotate("segment", x = x, xend = xend, y = y, yend = yend,
           colour = colour, linewidth = linewidth, linetype = linetype,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed"))
}
node <- function(x, y, w, h, fill = pal$light, colour = pal$ink, linewidth = 0.35) {
  annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
           fill = fill, colour = colour, linewidth = linewidth)
}

# Panel a: prior mechanism is solid; the present transcriptomic link is dashed.
p_a <- ggplot() +
  coord_cartesian(xlim = c(0, 14), ylim = c(0, 8.25), clip = "off") +
  node(1.05, 6.55, 1.7, 1.15) +
  node(3.55, 6.55, 2.7, 1.35, fill = pal$pale_blue, colour = pal$camp) +
  node(6.35, 6.55, 1.5, 1.15) +
  node(8.75, 6.55, 1.8, 1.15, fill = "#EDF4EE", colour = pal$camp) +
  node(11.55, 6.55, 2.25, 1.15) +
  arrow_segment(1.9, 2.2, 6.55, 6.55) +
  arrow_segment(4.9, 5.6, 6.55, 6.55) +
  arrow_segment(7.1, 7.85, 6.55, 6.55) +
  arrow_segment(9.7, 10.4, 6.55, 6.55) +
  node(3.6, 4.35, 3.35, 1.35, fill = "#F7F1E4", colour = pal$gold) +
  node(8.9, 4.35, 3.75, 1.35, fill = "#EAF3F4", colour = pal$teal) +
  arrow_segment(5.25, 7.0, 4.35, 4.35, colour = pal$teal, linetype = "22") +
  arrow_segment(7.05, 7.05, 4.95, 4.95, colour = pal$teal, linetype = "22") +
  annotate("text", x = 1.05, y = 6.55, label = "Nematode\ncues", size = 2.0, family = font_family) +
  annotate("text", x = 3.55, y = 6.78, label = "Reported GPCRs", size = 2.05, fontface = "bold", family = font_family) +
  annotate("text", x = 3.55, y = 6.28, label = "Gpr2 / Gpr3 / Gin3", size = 1.85, colour = pal$muted, family = font_family) +
  annotate("text", x = 6.35, y = 6.55, label = "Gpa2", size = 2.0, family = font_family) +
  annotate("text", x = 8.75, y = 6.55, label = "cAMP-PKA", size = 2.05, fontface = "bold", colour = pal$camp, family = font_family) +
  annotate("text", x = 11.55, y = 6.55, label = "Trap induction\nand development", size = 1.85, family = font_family) +
  annotate("text", x = 3.6, y = 4.65, label = "Prior evidence", size = 1.9, fontface = "bold", colour = pal$gold, family = font_family) +
  annotate("text", x = 3.6, y = 4.18, label = "GPCR loss changes Hog1\nphosphorylation dynamics", size = 1.72, colour = pal$ink, family = font_family) +
  annotate("text", x = 8.9, y = 4.65, label = "Present study", size = 1.9, fontface = "bold", colour = pal$teal, family = font_family) +
  annotate("text", x = 8.9, y = 4.18, label = "Pathway-mutant transcriptomes\n-> receptor mRNA response", size = 1.72, colour = pal$ink, family = font_family) +
  annotate("text", x = 6.8, y = 2.85, label = "Transcript-level feedback hypothesis", size = 1.9, fontface = "bold", colour = pal$teal, family = font_family) +
  annotate("text", x = 6.8, y = 2.38, label = "Three focal receptors: AOL_s00007g309,\nAOL_s00054g150 and AOL_s00176g42", size = 1.74, family = font_family) +
  arrow_segment(8.9, 6.8, 3.7, 3.12, colour = pal$teal, linetype = "22") +
  annotate("text", x = 0.15, y = 8.06, label = "Prior mechanism and the unanswered reverse direction", hjust = 0, size = 2.15, fontface = "bold", family = font_family) +
  annotate("text", x = 0.15, y = 7.72, label = "Solid arrows: reported mechanism; dashed arrows: transcriptomic association / hypothesis", hjust = 0, size = 1.78, colour = pal$muted, family = font_family) +
  theme_void(base_family = font_family) +
  theme(plot.margin = margin(2, 2, 2, 2))

# Panel b: a fixed eligibility chain, with the ranking used only after the final filter.
selection_counts[, x := seq_len(.N)]
selection_counts[, y := rev(seq_len(.N))]
selection_counts[, width := 1.35 + 0.072 * n]
selection_counts[, stage_short := c(
  "Complete GPCR models",
  "Prior 83-receptor bridge",
  "Prior 10 h induced",
  "Response in MAPK + cAMP-PKA",
  "Uncharacterized after prior work",
  "Focal receptors"
)]
p_b <- ggplot(selection_counts, aes(x = x, y = y)) +
  geom_tile(aes(x = 4.25 + width / 2, width = width, height = 0.62, fill = ifelse(x <= 4, "screen", "focus")), colour = "white", linewidth = 0.45) +
  geom_text(aes(x = 4.25 + width + 0.16, label = n), hjust = 0, size = 2.25, fontface = "bold", family = font_family, colour = pal$ink) +
  geom_text(aes(x = -2.35, label = stage_short), hjust = 0, size = 1.62, family = font_family, colour = pal$ink) +
  scale_fill_manual(values = c(screen = "#A8BAC9", focus = pal$camp), guide = "none") +
  scale_x_continuous(limits = c(-2.45, 13.7), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 6.75), breaks = NULL, expand = c(0, 0)) +
  annotate("text", x = 0.05, y = 6.58, label = "Fixed eligibility filters identify the focal set", hjust = 0, size = 2.1, fontface = "bold", family = font_family) +
  annotate("text", x = 0.05, y = 6.24, label = "No continuous score is used to define eligibility", hjust = 0, size = 1.72, colour = pal$muted, family = font_family) +
  annotate("text", x = 0.05, y = 0.52, label = "GIN2 / GIN6 / GIN9 are recovered as positive controls;\nwithin the final 13, effect size only sets follow-up order.", hjust = 0, size = 1.65, colour = pal$teal, family = font_family) +
  theme_void(base_family = font_family) +
  theme(plot.margin = margin(2, 2, 2, 2))

# Panel c: the quantitative hero panel, restricted to the three receptors.
heat_order <- c("AoHog1", "AoOpy2", "AoAcy", "AoPkaC1", "AoRab7", "AoVps35")
context[, regulator_label := factor(regulator_label, levels = heat_order)]
p_c <- ggplot(context, aes(x = regulator_label, y = display_id, fill = median_log2FC_mutant_vs_WT)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_tile(data = context[supported == TRUE], fill = NA, colour = pal$ink, linewidth = 0.72) +
  geom_text(aes(label = tile_label), size = 2.05, family = font_family, colour = pal$ink) +
  geom_vline(xintercept = c(2.5, 4.5), colour = "white", linewidth = 1.05) +
  scale_fill_gradient2(low = pal$down, mid = "#FAFAFA", high = pal$up, midpoint = 0,
                       limits = c(-3.5, 3.5), oob = squish, breaks = c(-3, -1.5, 0, 1.5, 3), name = "Median\nlog2FC") +
  scale_x_discrete(labels = c("AoHog1", "AoOpy2", "AoAcy", "AoPkaC1", "AoRab7", "AoVps35"), drop = FALSE) +
  labs(title = "Transcript responses are receptor- and module-specific",
       subtitle = "Matched-time mutant versus WT; * BH q < 0.05 in at least one matched contrast",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 6.6, base_family = font_family) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 5.65, angle = 38, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 5.75, colour = pal$ink), axis.title = element_blank(),
        plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_text(size = 5.85, colour = pal$muted, margin = margin(b = 3)),
        legend.position = "right", legend.title = element_text(size = 5.6), legend.text = element_text(size = 5.3),
        legend.key.height = unit(15, "mm"), plot.margin = margin(2, 2, 2, 2)) +
  annotate("text", x = 1.5, y = 3.68, label = "Osmosensing MAPK", size = 1.72, colour = pal$mapk, fontface = "bold", family = font_family) +
  annotate("text", x = 3.5, y = 3.68, label = "cAMP-PKA", size = 1.72, colour = pal$camp, fontface = "bold", family = font_family) +
  annotate("text", x = 5.5, y = 3.68, label = "Vacuolar trafficking", size = 1.72, colour = pal$vac, fontface = "bold", family = font_family)

# Panel d: exposure-linked prior evidence, with heterogeneous induction made explicit.
p_d_ratio <- ggplot(exposure, aes(x = prior_log2ratio, y = display_id, colour = direction)) +
  geom_vline(xintercept = 0, colour = pal$grey, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = prior_log2ratio, yend = display_id), colour = pal$grey, linewidth = 0.65) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%+.2f", prior_log2ratio)), hjust = ifelse(exposure$prior_log2ratio >= 0, -0.28, 1.2), size = 1.85, colour = pal$ink, family = font_family) +
  scale_colour_manual(values = c(Up = pal$up, Down = pal$down), guide = "none") +
  scale_x_continuous(limits = c(-0.9, 7.2), breaks = c(0, 2, 4, 6), expand = c(0, 0)) +
  labs(title = "Prior TWF154 exposure snapshot", subtitle = "Mean TPM, 10 h versus 0 h", x = "log2 ratio", y = NULL) +
  theme_classic(base_size = 6.3, base_family = font_family) +
  theme(plot.title = element_text(size = 6.8, face = "bold", margin = margin(b = 1)), plot.subtitle = element_text(size = 5.45, colour = pal$muted, margin = margin(b = 2)), axis.text.y = element_text(size = 5.35), axis.text.x = element_text(size = 5.25), axis.title.x = element_text(size = 5.65), plot.margin = margin(2, 3, 1, 2))

p_d_support <- ggplot(exposure, aes(x = cross_geo_support, y = display_id)) +
  geom_col(width = 0.55, fill = pal$camp) +
  geom_text(aes(label = sprintf("%d/%d", cross_geo_support, cross_geo_total)), hjust = -0.15, size = 1.85, family = font_family, colour = pal$ink) +
  scale_x_continuous(limits = c(0, 7.3), breaks = c(0, 2, 4, 6), expand = c(0, 0)) +
  labs(title = "ATCC cross-GEO induction support", subtitle = "Independent series with supported WT induction", x = "Series", y = NULL) +
  theme_classic(base_size = 6.3, base_family = font_family) +
  theme(plot.title = element_text(size = 6.8, face = "bold", margin = margin(b = 1)), plot.subtitle = element_text(size = 5.45, colour = pal$muted, margin = margin(b = 2)), axis.text.y = element_text(size = 5.35), axis.text.x = element_text(size = 5.25), axis.title.x = element_text(size = 5.65), plot.margin = margin(2, 3, 1, 2))

p_d_stack <- p_d_ratio / p_d_support + plot_layout(heights = c(1, 0.95), tag_level = "keep")
p_d <- wrap_elements(full = p_d_stack)

figure <- (p_a + p_b + plot_layout(widths = c(1.33, 0.87))) /
  (p_c + p_d + plot_layout(widths = c(1.38, 0.82))) +
  plot_layout(heights = c(0.9, 1.35)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(size = 8.2, face = "bold", family = font_family), plot.tag.position = c(0, 1), plot.margin = margin(2, 2, 2, 2)))

base_out <- file.path(figure_dir, "Fig_opening_three_receptors_context")
width_mm <- 183
height_mm <- 145
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(paste0(base_out, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
grDevices::cairo_pdf(paste0(base_out, ".pdf"), width = width_in, height = height_in, family = font_family)
print(figure)
dev.off()
ragg::agg_tiff(paste0(base_out, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(figure)
dev.off()
ragg::agg_png(paste0(base_out, ".png"), width = width_in, height = height_in, units = "in", res = 300)
print(figure)
dev.off()

cat("Figure:", paste0(base_out, ".{svg,pdf,tiff,png}"), "\n")
cat("Focal context rows:", nrow(context), "\n")
cat("Selection stages:", paste(selection_counts$n, collapse = " -> "), "\n")
cat("Prior exposure ratios:", paste(sprintf("%s=%+.2f", exposure$receptor, exposure$prior_log2ratio), collapse = "; "), "\n")
