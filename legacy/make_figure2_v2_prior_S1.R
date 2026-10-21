#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "Usage: Rscript make_figure2_v2.R",
    "<WT_atlas.tsv> <evidence_by_gse.tsv> <tree_metadata.tsv>",
    "<selection.tsv> <focal_response.tsv> <design_manifest.tsv> <output_dir>"
  ))
}

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

atlas_path <- normalizePath(args[[1]], mustWork = TRUE)
evidence_path <- normalizePath(args[[2]], mustWork = TRUE)
tree_path <- normalizePath(args[[3]], mustWork = TRUE)
selection_path <- normalizePath(args[[4]], mustWork = TRUE)
response_path <- normalizePath(args[[5]], mustWork = TRUE)
design_path <- normalizePath(args[[6]], mustWork = TRUE)
output_dir <- args[[7]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             encoding = "UTF-8")
}

atlas_raw <- read_tsv(atlas_path)
evidence_raw <- read_tsv(evidence_path)
tree <- read_tsv(tree_path)
selection <- read_tsv(selection_path)
response <- read_tsv(response_path)
design <- read_tsv(design_path)

stopifnot(all(c("candidate_id", "gse", "log2fc") %in% names(atlas_raw)))
stopifnot(all(c("candidate_id", "gse", "wt_induction_supported") %in%
              names(evidence_raw)))
stopifnot("gpcr" %in% names(tree))
stopifnot(all(c("candidate_id", "focal", "Kuo_gene_name",
                "max_abs_l2_MAPK", "max_abs_l2_cAMP_PKA", "sum") %in%
              names(selection)))
stopifnot(all(c("receptor", "AOL_id", "regulator_background",
                "pathway_module", "median_log2FC_mutant_vs_WT",
                "n_supported_BHq005") %in% names(response)))
stopifnot(all(c("matrix_id", "gse", "unit", "genotypes", "time_points",
                "treatments", "matrix_status") %in% names(design)))

is_true <- function(x) tolower(as.character(x)) == "true"
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

inventory <- unique(tree$gpcr)
stopifnot(length(inventory) == 89L, !anyNA(inventory), !anyDuplicated(inventory))

gse_order <- c(
  "GSE210886", "GSE212735", "GSE213542", "GSE215079", "GSE245014",
  "GSE253759", "GSE301423", "GSE301452", "GSE303416", "GSE306156",
  "GSE322548"
)
stopifnot(all(gse_order %in% unique(atlas_raw$gse)))

gse_display <- c(
  GSE210886 = "GSE210886\nAcy/PKA",
  GSE212735 = "GSE212735\nAtg11/N-starv.",
  GSE213542 = "GSE213542\nbZIP129",
  GSE215079 = "GSE215079\nAge2",
  GSE245014 = "GSE245014\nRab7/Vps35",
  GSE253759 = "GSE253759\nAge2*",
  GSE301423 = "GSE301423\nHam4/Ham8",
  GSE301452 = "GSE301452\nChk1",
  GSE303416 = "GSE303416\nHal4b†",
  GSE306156 = "GSE306156\nOpy2",
  GSE322548 = "GSE322548\nHog1"
)

selection$focal_flag <- is_true(selection$focal)
selection$known_flag <- nzchar(trimws(selection$Kuo_gene_name))
focal_ids <- selection$candidate_id[selection$focal_flag]
known_ids <- selection$candidate_id[selection$known_flag]
uncharacterized <- selection[!selection$known_flag, , drop = FALSE]
uncharacterized <- uncharacterized[order(-uncharacterized$sum), , drop = FALSE]

stopifnot(nrow(selection) == 16L)
stopifnot(length(focal_ids) == 3L)
stopifnot(length(known_ids) == 3L)
stopifnot(nrow(uncharacterized) == 13L)
stopifnot(all(selection$candidate_id %in% inventory))
stopifnot(identical(
  uncharacterized$candidate_id[seq_len(3L)],
  c("AOL_s00007g309.1", "AOL_s00054g150.1", "AOL_s00176g42.1")
))

atlas_ids <- unique(atlas_raw$candidate_id)
inventory_mapped <- intersect(inventory, atlas_ids)
inventory_unmapped <- setdiff(inventory, atlas_ids)
atlas_outside_inventory <- setdiff(atlas_ids, inventory)
stopifnot(length(atlas_ids) == 90L)
stopifnot(length(inventory_mapped) == 87L)
stopifnot(length(inventory_unmapped) == 2L)
stopifnot(length(atlas_outside_inventory) == 3L)

rank_map <- tapply(atlas_raw$priority_rank, atlas_raw$candidate_id,
                   function(z) z[which(!is.na(z))[1]])
other_ids <- setdiff(inventory, selection$candidate_id)
other_rank <- rank_map[other_ids]
other_ids <- other_ids[order(is.na(other_rank), other_rank, other_ids)]
selected_other_ids <- uncharacterized$candidate_id[
  !uncharacterized$candidate_id %in% focal_ids
]
known_ids <- selection$candidate_id[selection$known_flag][
  order(-selection$sum[selection$known_flag])
]
row_order <- c(
  c("AOL_s00007g309.1", "AOL_s00054g150.1", "AOL_s00176g42.1"),
  selected_other_ids,
  known_ids,
  other_ids
)
stopifnot(length(row_order) == 89L, !anyDuplicated(row_order),
          setequal(row_order, inventory))

grid <- expand.grid(candidate_id = row_order, gse = gse_order,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
atlas_sub <- atlas_raw[
  atlas_raw$candidate_id %in% inventory & atlas_raw$gse %in% gse_order,
  c("candidate_id", "gse", "log2fc", "priority_rank", "priority_tier")
]
grid <- merge(grid, atlas_sub, by = c("candidate_id", "gse"), all.x = TRUE,
              sort = FALSE)
evidence_sub <- evidence_raw[
  evidence_raw$candidate_id %in% inventory & evidence_raw$gse %in% gse_order,
  c("candidate_id", "gse", "wt_induction_records",
    "wt_induction_supported", "wt_induction_direction")
]
grid <- merge(grid, evidence_sub, by = c("candidate_id", "gse"), all.x = TRUE,
              sort = FALSE)
grid$row_index <- match(grid$candidate_id, row_order)
grid$column_index <- match(grid$gse, gse_order)
grid <- grid[order(grid$row_index, grid$column_index), , drop = FALSE]
grid$status <- "Other GPCR"
grid$status[grid$candidate_id %in% selected_other_ids] <-
  "Filter-positive uncharacterized"
grid$status[grid$candidate_id %in% known_ids] <- "Previously named/tested"
grid$status[grid$candidate_id %in% focal_ids] <- "Focal GPCR"
grid$within_figure1_inventory <- TRUE
grid$mapped_in_cross_gse_atlas <- grid$candidate_id %in% atlas_ids
grid$independent_series <- grid$gse != "GSE253759"
grid$provisional_mapping <- grid$gse == "GSE303416"

stopifnot(nrow(grid) == 89L * 11L)
stopifnot(sum(is.na(grid$log2fc)) == 2L * 11L)

response_focus_names <- c(
  "AOL_s00007g309", "AOL_s00054g150", "AOL_s00176g42"
)
response_focus <- response[response$receptor %in% response_focus_names, ,
                           drop = FALSE]
display_background <- c(
  "Kuo et al. WT 0→10 hpe" = "0–10 hpe",
  "ΔAoAcy" = "ΔAoAcy",
  "ΔAoPkaC1" = "ΔAoPkaC1",
  "ΔAoOpy2" = "ΔAoOpy2",
  "ΔAoHog1" = "ΔAoHog1",
  "ΔAoRab7" = "ΔAoRab7",
  "ΔAoVps35" = "ΔAoVps35"
)
background_order <- names(display_background)
response_focus <- response_focus[
  match(as.vector(outer(response_focus_names, background_order,
                        paste, sep = "||")),
        paste(response_focus$receptor, response_focus$regulator_background,
              sep = "||")),
  , drop = FALSE
]
stopifnot(!anyNA(response_focus$receptor), nrow(response_focus) == 21L)
response_focus$supported <- !is.na(response_focus$n_supported_BHq005) &
  response_focus$n_supported_BHq005 > 0
response_focus$display_background <- unname(
  display_background[response_focus$regulator_background]
)

design_wt <- design[design$gse %in% gse_order & design$unit != "count", ,
                    drop = FALSE]
design_wt$included_in_panel_a <- TRUE
design_day5 <- design[design$gse == "GSE322740", , drop = FALSE]
design_day5$included_in_panel_a <- FALSE
design_export <- rbind(
  design_wt[, intersect(names(design_wt), c(
    "matrix_id", "gse", "unit", "genotypes", "time_points", "treatments",
    "mapping_confidence", "matrix_status", "notes", "included_in_panel_a"
  ))],
  design_day5[, intersect(names(design_day5), c(
    "matrix_id", "gse", "unit", "genotypes", "time_points", "treatments",
    "mapping_confidence", "matrix_status", "notes", "included_in_panel_a"
  ))]
)

write.table(grid, file.path(output_dir, "Figure2_SourceData_PanelA.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(selection, file.path(output_dir, "Figure2_SourceData_PanelB.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(uncharacterized,
            file.path(output_dir, "Figure2_SourceData_PanelC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(response_focus,
            file.path(output_dir, "Figure2_SourceData_PanelD.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(design_export,
            file.path(output_dir, "Figure2_SourceData_StudyDesign.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

mapping_audit <- data.frame(
  category = c(
    "Figure1_inventory", "Cross_GEO_atlas", "Intersection",
    "Figure1_only_unmapped", "Atlas_only_excluded"
  ),
  n = c(89L, 90L, 87L, length(inventory_unmapped),
        length(atlas_outside_inventory)),
  ids = c(
    "", "", "", paste(inventory_unmapped, collapse = ";"),
    paste(atlas_outside_inventory, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.table(mapping_audit,
            file.path(output_dir, "Figure2_candidate_mapping_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

pal <- c(
  ink = "#232323",
  neutral = "#BFC4C9",
  neutral_light = "#ECEEEF",
  selected = "#74B8B3",
  selected_dark = "#2F7F7A",
  focal = "#D45A4A",
  focal_light = "#F4C9C2",
  known = "#6F617B",
  mapk = "#4E79A7",
  pka = "#E39C45",
  negative = "#326AA5",
  zero = "#F8F8F5",
  positive = "#C84B4B",
  missing = "#D9DDE0"
)

heat_cols <- colorRampPalette(c(pal[["negative"]], pal[["zero"]],
                               pal[["positive"]]))(201)
zlim <- c(-5, 5)
color_index <- function(z) {
  z <- pmin(zlim[2], pmax(zlim[1], z))
  1L + round((z - zlim[1]) / diff(zlim) * 200)
}

panel_header <- function(letter, title) {
  mtext(paste0(letter, "   ", title), side = 3, line = 0.95, adj = 0,
        font = 2, cex = 0.90, col = pal[["ink"]])
}

draw_box <- function(xleft, ybottom, xright, ytop, label,
                     fill = "white", border = pal[["neutral"]],
                     text_col = pal[["ink"]], font = 1, cex = 0.72) {
  rect(xleft, ybottom, xright, ytop, col = fill, border = border, lwd = 0.8)
  text((xleft + xright) / 2, (ybottom + ytop) / 2, label,
       col = text_col, font = font, cex = cex)
}

plot_panel_a <- function() {
  nr <- length(row_order)
  nc <- length(gse_order)
  par(mar = c(6.4, 1.2, 2.4, 7.2), xpd = NA)
  plot.new()
  plot.window(xlim = c(0.25, nc + 2.85), ylim = c(-9.0, nr + 1.2),
              xaxs = "i", yaxs = "i")
  panel_header("a", "Cross-GEO WT-response atlas of 89 GPCRs")

  for (i in seq_len(nrow(grid))) {
    x <- grid$column_index[i]
    y <- nr - grid$row_index[i] + 1
    val <- grid$log2fc[i]
    fill <- if (is.na(val)) pal[["missing"]] else heat_cols[color_index(val)]
    rect(x - 0.48, y - 0.46, x + 0.48, y + 0.46,
         col = fill, border = "white", lwd = 0.16)
    if (!is.na(grid$wt_induction_supported[i]) &&
        grid$wt_induction_supported[i] > 0) {
      points(x, y, pch = 16, cex = 0.19,
             col = if (!is.na(val) && abs(val) >= 3.2) "white" else pal[["ink"]])
    }
  }

  status_by_row <- unique(grid[, c("candidate_id", "row_index", "status",
                                   "mapped_in_cross_gse_atlas")])
  status_by_row <- status_by_row[order(status_by_row$row_index), ]
  status_cols <- c(
    "Other GPCR" = pal[["neutral_light"]],
    "Filter-positive uncharacterized" = pal[["selected"]],
    "Previously named/tested" = pal[["known"]],
    "Focal GPCR" = pal[["focal"]]
  )
  for (i in seq_len(nrow(status_by_row))) {
    y <- nr - status_by_row$row_index[i] + 1
    rect(0.38, y - 0.46, 0.72, y + 0.46,
         col = status_cols[status_by_row$status[i]], border = NA)
  }

  selected_cut_y <- nr - 16 + 0.5
  segments(0.32, selected_cut_y, nc + 0.49, selected_cut_y,
           lwd = 0.85, col = pal[["ink"]])
  text(0.26, nr - 7.5, "16 filter-positive", srt = 90,
       cex = 0.63, font = 2, col = pal[["ink"]])

  focal_label_y <- nr - c(0.0, 2.8, 5.6)
  for (j in seq_along(focal_ids)) {
    id <- focal_ids[j]
    row_i <- match(id, row_order)
    y <- nr - row_i + 1
    segments(nc + 0.49, y, nc + 0.70, y, col = pal[["focal"]], lwd = 0.65)
    segments(nc + 0.70, y, nc + 0.90, focal_label_y[j],
             col = pal[["focal"]], lwd = 0.65)
    text(nc + 0.96, focal_label_y[j], strip_version(id), pos = 4,
         cex = 0.54, font = 2, col = pal[["focal"]])
  }

  known_rows <- match(known_ids, row_order)
  known_y <- nr - known_rows + 1
  segments(nc + 0.54, min(known_y), nc + 0.54, max(known_y),
           col = pal[["known"]], lwd = 0.7)
  segments(nc + 0.49, min(known_y), nc + 0.66, min(known_y),
           col = pal[["known"]], lwd = 0.7)
  segments(nc + 0.49, max(known_y), nc + 0.66, max(known_y),
           col = pal[["known"]], lwd = 0.7)
  text(nc + 0.76, mean(known_y), "GIN2 / GIN6 / GIN9", pos = 4,
       cex = 0.54, font = 2, col = pal[["known"]])

  unmapped_rows <- match(inventory_unmapped, row_order)
  unmapped_y <- nr - unmapped_rows + 1
  segments(nc + 0.54, min(unmapped_y), nc + 0.54, max(unmapped_y),
           col = "#6D7378", lwd = 0.65)
  segments(nc + 0.49, min(unmapped_y), nc + 0.66, min(unmapped_y),
           col = "#6D7378", lwd = 0.65)
  segments(nc + 0.49, max(unmapped_y), nc + 0.66, max(unmapped_y),
           col = "#6D7378", lwd = 0.65)
  text(nc + 0.76, mean(unmapped_y), "2 Figure 1 entries not mapped",
       pos = 4, cex = 0.52, col = "#6D7378")

  text(seq_len(nc), -0.55, unname(gse_display[gse_order]), srt = 55,
       adj = 1, cex = 0.53, col = pal[["ink"]])

  key_y <- seq(35, 59, length.out = 51)
  for (i in seq_len(length(key_y) - 1L)) {
    idx <- round((i - 1) / (length(key_y) - 2) * 200) + 1L
    rect(nc + 1.02, key_y[i], nc + 1.24, key_y[i + 1],
         col = heat_cols[idx], border = NA)
  }
  text(nc + 1.32, c(35, 47, 59), c("−5", "0", "+5"),
       pos = 4, cex = 0.52)
  text(nc + 1.13, 62.0, "WT log2FC", cex = 0.54, font = 2)
  points(nc + 1.12, 31.8, pch = 16, cex = 0.20, col = pal[["ink"]])
  text(nc + 1.34, 31.8, "q < 0.05", pos = 4, cex = 0.52)

  legend_y <- c(27.6, 24.5, 21.4, 18.3, 15.2)
  legend_labels <- c("Focal", "Other selected", "Named/tested", "Other",
                     "Not mapped")
  legend_fills <- c(pal[["focal"]], pal[["selected"]], pal[["known"]],
                    pal[["neutral_light"]], pal[["missing"]])
  for (i in seq_along(legend_y)) {
    rect(nc + 1.01, legend_y[i] - 0.38, nc + 1.24,
         legend_y[i] + 0.38, col = legend_fills[i], border = NA)
    text(nc + 1.34, legend_y[i], legend_labels[i], pos = 4, cex = 0.50)
  }
}

plot_panel_b <- function() {
  par(mar = c(0.7, 0.5, 2.4, 0.5), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  panel_header("b", "Predefined prioritization")

  draw_box(0.10, 0.80, 0.90, 0.96,
           "89 GPCRs from\nATCC reference genome",
           fill = "#F4F5F6", border = pal[["neutral"]], cex = 0.72)
  arrows(0.50, 0.79, 0.50, 0.72, length = 0.045, lwd = 0.75)
  draw_box(0.13, 0.57, 0.87, 0.71,
           "16 met four binary\ntranscriptome criteria",
           fill = "#DDECEB", border = pal[["selected_dark"]],
           font = 2, cex = 0.72)
  arrows(0.50, 0.56, 0.50, 0.49, length = 0.045, lwd = 0.75)
  draw_box(0.09, 0.33, 0.91, 0.48,
           "Remove named/tested\nGIN2, GIN6 and GIN9",
           fill = "#F1EEF3", border = pal[["known"]], cex = 0.70)
  arrows(0.50, 0.32, 0.50, 0.25, length = 0.045, lwd = 0.75)
  draw_box(0.13, 0.12, 0.87, 0.24,
           "13 uncharacterized\ncandidates ranked",
           fill = "white", border = pal[["selected_dark"]], cex = 0.70)
  arrows(0.50, 0.11, 0.50, 0.065, length = 0.040, lwd = 0.75)
  text(0.50, 0.015, "Advance top three", font = 2, cex = 0.74,
       col = pal[["focal"]])
}

plot_panel_c <- function() {
  par(mar = c(3.8, 7.4, 2.4, 0.9), xpd = NA)
  score_mat <- rbind(uncharacterized$max_abs_l2_MAPK,
                     uncharacterized$max_abs_l2_cAMP_PKA)
  display_order <- rev(seq_len(nrow(uncharacterized)))
  score_mat <- score_mat[, display_order, drop = FALSE]
  dat <- uncharacterized[display_order, , drop = FALSE]
  mids <- barplot(score_mat, horiz = TRUE, beside = FALSE,
                  col = c(pal[["mapk"]], pal[["pka"]]), border = NA,
                  space = 0.22, xlim = c(0, 7.0), axes = FALSE,
                  names.arg = rep("", nrow(dat)))
  axis(1, at = 0:7, labels = 0:7, cex.axis = 0.73,
       col = pal[["ink"]], col.axis = pal[["ink"]], lwd = 0,
       lwd.ticks = 0.55)
  abline(v = 1:6, col = "#ECEDEF", lwd = 0.5)
  panel_header("c", "Ranking of 13 candidates")
  mtext("Combined maximum absolute log2FC", side = 1,
        line = 2.25, cex = 0.73)

  short_ids <- strip_version(dat$candidate_id)
  focal_here <- dat$candidate_id %in% focal_ids
  usr <- par("usr")
  text(usr[1] - 0.10, mids, short_ids, pos = 2, xpd = NA,
       cex = 0.62, font = ifelse(focal_here, 2, 1),
       col = ifelse(focal_here, pal[["focal"]], pal[["ink"]]))
  text(dat$sum + 0.07, mids, sprintf("%.2f", dat$sum), pos = 4,
       cex = 0.58, col = pal[["ink"]], xpd = NA)
  for (i in which(focal_here)) {
    rect(0, mids[i] - 0.43, dat$sum[i], mids[i] + 0.43,
         border = pal[["focal"]], lwd = 0.9)
  }
  legend("bottomright", inset = c(0.00, 0.00), bty = "n", cex = 0.61,
         fill = c(pal[["mapk"]], pal[["pka"]]), border = NA,
         legend = c("Osmosensing MAPK", "cAMP–PKA"))
}

plot_panel_d <- function() {
  par(mar = c(5.1, 6.1, 2.4, 0.5), xpd = NA)
  nr <- length(response_focus_names)
  nc <- length(background_order)
  value_mat <- matrix(NA_real_, nrow = nr, ncol = nc,
                      dimnames = list(response_focus_names,
                                      unname(display_background)))
  support_mat <- matrix(FALSE, nrow = nr, ncol = nc,
                        dimnames = dimnames(value_mat))
  for (i in seq_len(nrow(response_focus))) {
    r <- match(response_focus$receptor[i], response_focus_names)
    c <- match(response_focus$regulator_background[i], background_order)
    value_mat[r, c] <- response_focus$median_log2FC_mutant_vs_WT[i]
    support_mat[r, c] <- isTRUE(response_focus$supported[i])
  }

  plot.new()
  plot.window(xlim = c(0.5, nc + 0.5), ylim = c(-1.15, nr + 0.95),
              xaxs = "i", yaxs = "i")
  panel_header("d", "Focal response profiles")
  y_positions <- rev(seq_len(nr))
  for (r in seq_len(nr)) {
    y <- y_positions[r]
    for (c in seq_len(nc)) {
      val <- value_mat[r, c]
      rect(c - 0.48, y - 0.48, c + 0.48, y + 0.48,
           col = heat_cols[color_index(val)], border = "white", lwd = 0.7)
      text(c, y, sprintf("%.1f", val), cex = 0.62,
           col = if (abs(val) >= 2.8) "white" else pal[["ink"]],
           font = 2)
      if (support_mat[r, c]) {
        points(c + 0.32, y + 0.31, pch = 16, cex = 0.38,
               col = if (abs(val) >= 2.8) "white" else pal[["ink"]])
      }
    }
  }
  text(0.38, y_positions, response_focus_names, pos = 2, cex = 0.61,
       font = 2, xpd = NA)
  text(seq_len(nc), 0.37, unname(display_background), srt = 48,
       adj = 1, cex = 0.60, xpd = NA)

  key_x <- seq(1.3, 5.0, length.out = 41)
  for (i in seq_len(length(key_x) - 1L)) {
    idx <- round((i - 1) / (length(key_x) - 2) * 200) + 1L
    rect(key_x[i], -0.88, key_x[i + 1], -0.70,
         col = heat_cols[idx], border = NA)
  }
  text(c(1.3, 3.15, 5.0), -0.99, c("−5", "0", "+5"), cex = 0.56)
  text(3.15, -0.50, "Signed log2 expression change", cex = 0.59,
       font = 2)
  points(5.55, -0.80, pch = 16, cex = 0.38, col = pal[["ink"]])
  text(5.72, -0.80, "Supported", pos = 4, cex = 0.57)
}

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(family = "Arial", fg = pal[["ink"]], col = pal[["ink"]],
      col.axis = pal[["ink"]], col.lab = pal[["ink"]],
      lend = "round", ljoin = "round", ps = 10)
  layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE),
         widths = c(0.78, 1.42, 1.10), heights = c(1.18, 0.82))
  plot_panel_a()
  plot_panel_b()
  plot_panel_c()
  plot_panel_d()
}

width_mm <- 183
height_mm <- 205
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
base_name <- file.path(output_dir, "Figure2_revised")

grDevices::cairo_pdf(paste0(base_name, ".pdf"), width = width_in,
                     height = height_in, family = "Arial", onefile = TRUE)
draw_figure()
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(base_name, ".svg"), width = width_in,
                   height = height_in,
                   system_fonts = list(sans = "Arial"))
} else {
  grDevices::svg(paste0(base_name, ".svg"), width = width_in,
                 height = height_in, family = "Arial", onefile = TRUE)
}
draw_figure()
dev.off()

grDevices::png(paste0(base_name, ".png"), width = width_in,
               height = height_in, units = "in", res = 600,
               type = "cairo", family = "Arial")
draw_figure()
dev.off()

grDevices::png(file.path(output_dir, "Figure2_revised_preview.png"),
               width = width_in, height = height_in, units = "in", res = 300,
               type = "cairo", family = "Arial")
draw_figure()
dev.off()

qa_dir <- file.path(output_dir, "qa_panels")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
save_qa_panel <- function(plot_fun, name, width, height) {
  grDevices::png(file.path(qa_dir, paste0(name, ".png")),
                 width = width, height = height, units = "in", res = 300,
                 type = "cairo", family = "Arial")
  par(family = "Arial", fg = pal[["ink"]], col = pal[["ink"]],
      col.axis = pal[["ink"]], col.lab = pal[["ink"]],
      lend = "round", ljoin = "round", ps = 10)
  plot_fun()
  dev.off()
}
save_qa_panel(plot_panel_a, "panel_a", width_in, 4.6)
save_qa_panel(plot_panel_b, "panel_b", 2.2, 3.6)
save_qa_panel(plot_panel_c, "panel_c", 3.5, 3.6)
save_qa_panel(plot_panel_d, "panel_d", 3.2, 3.6)

audit_lines <- c(
  "Figure 2 revised data audit",
  paste("Figure 1 inventory:", length(inventory)),
  paste("Cross-GEO atlas candidates:", length(atlas_ids)),
  paste("Intersection retained in panel A:", length(inventory_mapped)),
  paste("Figure 1 candidates not mapped in atlas:",
        paste(inventory_unmapped, collapse = ", ")),
  paste("Atlas candidates excluded because absent from Figure 1 inventory:",
        paste(atlas_outside_inventory, collapse = ", ")),
  paste("Panel A cells:", nrow(grid)),
  paste("Panel A missing cells:", sum(is.na(grid$log2fc))),
  paste("Panel A GEO series:", length(gse_order)),
  "GSE253759 is displayed but flagged as non-independent from GSE215079.",
  "GSE303416 is displayed and flagged as having provisional/order-based mapping.",
  "GSE322740 is excluded from panel A because the Day 5 design has no WT time-induction contrast.",
  paste("Panel B filter-positive candidates:", nrow(selection)),
  paste("Panel C uncharacterized candidates:", nrow(uncharacterized)),
  paste("Panel D signed focal-response cells:", nrow(response_focus)),
  paste("Focal candidates:", paste(strip_version(focal_ids), collapse = ", ")),
  paste("Source MD5 [atlas]:", unname(tools::md5sum(atlas_path))),
  paste("Source MD5 [evidence]:", unname(tools::md5sum(evidence_path))),
  paste("Source MD5 [tree]:", unname(tools::md5sum(tree_path))),
  paste("Source MD5 [selection]:", unname(tools::md5sum(selection_path))),
  paste("Source MD5 [response]:", unname(tools::md5sum(response_path))),
  paste("Source MD5 [design]:", unname(tools::md5sum(design_path)))
)
writeLines(audit_lines, file.path(output_dir, "Figure2_revised_data_audit.txt"))
