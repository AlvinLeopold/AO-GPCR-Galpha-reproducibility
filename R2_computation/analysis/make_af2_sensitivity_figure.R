#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript make_af2_sensitivity_figure.R <analysis_dir> <output_dir>")
analysis_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.delim(file.path(analysis_dir, "af2_all_model_scores_clean.tsv"), check.names = FALSE)
rob <- read.delim(file.path(analysis_dir, "af2_receptor_rank_robustness.tsv"), check.names = FALSE)
focal <- c(AoGPR200 = "AOL_s00007g309.1", AoGPRPR8 = "AOL_s00176g42.1")
col <- c(AOGpa1 = "#3A668E", AOGpa2 = "#C76A4D", AOGpa3 = "#579F91")

draw <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(1:3, nrow = 1), widths = c(1.05, 1.05, 1.2))
  par(mar = c(4.5, 4.5, 2.8, 0.8), family = "Arial", las = 1)
  for (i in seq_along(focal)) {
    s <- scores[scores$gpcr == focal[[i]], ]
    plot(NA, xlim = c(0.5, 3.5), ylim = c(0, 0.85), xlab = "Partner",
         ylab = if (i == 1L) "ipTM" else "", xaxt = "n", main = names(focal)[i])
    axis(1, at = 1:3, labels = c("AoGpa1", "AoGpa2", "AoGpa3"))
    abline(h = seq(0, 0.8, 0.2), col = "#E8ECF0", lwd = 0.8)
    for (k in 1:3) {
      part <- names(col)[k]
      v <- s[s$galpha == part, ]
      points(rep(k, nrow(v)) + c(-0.10, 0, 0.10)[v$model], v$iptm,
             pch = 16, cex = 1.1, col = col[[part]])
      segments(k - 0.17, median(v$iptm), k + 0.17, median(v$iptm),
               lwd = 2.5, col = "#202C35")
    }
    mtext("3 model weights, one seed", side = 3, line = 0.1, cex = 0.72, col = "#525E68")
  }
  par(mar = c(4.5, 4.8, 2.8, 1), family = "Arial", las = 1)
  changed <- !as.logical(rob$max_mean_agree)
  plot(rob$max_margin, rob$mean_margin, xlim = c(0, max(rob$max_margin) * 1.04),
       ylim = c(0, max(rob$mean_margin) * 1.04), pch = 16, cex = 0.62,
       col = ifelse(changed, "#C76A4D", "#59798D66"),
       xlab = "Top-partner margin from max ipTM",
       ylab = "Top-partner margin from mean ipTM", main = "89-receptor ranking")
  abline(v = c(0.02, 0.05), lty = 3, col = "#B3BDC4")
  abline(0, 1, lty = 2, col = "#A5AFB7")
  legend("topleft", legend = c("Same top partner", "Top partner changes"),
         col = c("#59798D", "#C76A4D"), pch = 16, cex = 0.72, bty = "n")
}

pdf_path <- file.path(output_dir, "Supplementary_Figure_AF2_model_variability.pdf")
png_path <- file.path(output_dir, "Supplementary_Figure_AF2_model_variability.png")
grDevices::cairo_pdf(pdf_path, width = 11.5, height = 4.0, family = "Arial")
draw(); dev.off()
grDevices::png(png_path, width = 11.5, height = 4.0, units = "in", res = 300,
               type = "cairo", family = "Arial")
draw(); dev.off()
cat(pdf_path, "\n", png_path, "\n")
