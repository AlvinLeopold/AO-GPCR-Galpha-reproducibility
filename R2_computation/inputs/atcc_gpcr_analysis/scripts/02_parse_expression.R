#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})
source(file.path(getwd(), "scripts", "common.R"), local = TRUE)

step <- "02_parse_expression"
candidates_path <- file.path(PROJ, "metadata", "candidates.tsv")
cfg_path <- file.path(PROJ, "metadata", "source_file_manifest.tsv")
catalog_sample_path <- file.path(REPO, "geo_ao_transcriptomes", "catalog", "sample_metadata.tsv")
start_log(step, c(candidates_path, cfg_path, catalog_sample_path))

candidates <- fread(candidates_path, sep = "\t", encoding = "UTF-8")
cfg <- fread(cfg_path, sep = "\t", encoding = "UTF-8")
catalog_samples <- fread(catalog_sample_path, sep = "\t", encoding = "UTF-8")

## Every expression column is specified explicitly. This keeps sample identity
## independent of column order and records the few GEO naming corrections.
make_spec <- function(matrix_id, gse, unit, file, cols, titles, notes = "") {
  if (length(cols) != length(titles)) stop("Column/title length mismatch for ", matrix_id)
  list(matrix_id = matrix_id, gse = gse, unit = unit, file = file,
       cols = cols, titles = titles, notes = notes)
}

cfg_file <- function(source_id) file.path(PROJ, cfg$file[match(source_id, cfg$source_id)])

specs <- list(
  make_spec("GSE210886_Acy", "GSE210886", "TPM", cfg_file("GSE210886_Acy"),
            c(paste0("AO_0_", 1:3), paste0("AO_12_", 1:3), paste0("Acy_0_", 1:3), paste0("Acy_12_", 1:3)),
            c(paste0("AO_0_", 1:3), paste0("AO_12_", 1:3), paste0("Acy_0_", 1:3), paste0("Acy_12_", 1:3)),
            "GEO sample titles are identical to matrix columns; AO is WT and Acy is delta-AoAcy."),
  make_spec("GSE210886_PKA", "GSE210886", "TPM", cfg_file("GSE210886_PKA"),
            c(paste0("A_0_", 1:3), paste0("A_12_", 1:3), paste0("C_0_", 1:3), paste0("C_12_", 1:3)),
            c(paste0("A_0_", 1:3), paste0("A_12_", 1:3), paste0("C_0_", 1:3), paste0("C_12_", 1:3)),
            "The A/C block is a distinct WT/delta-AoPkaC1 library set and is not deduplicated with AO/Acy."),
  make_spec("GSE212735", "GSE212735", "TPM", cfg_file("GSE212735"),
            c(paste0("A11_0_", 1:3), paste0("A11_N_", 1:3), paste0("AO_0_", 1:3), paste0("WT_N_", 1:3)),
            c(paste0("A11_0_", 1:3), paste0("A11_N_", 1:3), paste0("A_0_", 1:3), paste0("WT_N_", 1:3)),
            "The GEO sample titles call the WT baseline A_0; the matrix uses AO_0."),
  make_spec("GSE213542", "GSE213542", "FPKM", cfg_file("GSE213542"),
            c(paste0("WT_0_", 1:3), paste0("WT_1_", 1:3), paste0("Ao129_0_", 1:3), paste0("Ao129_1_", 1:3)),
            c(paste0("WT-0-", 1:3), paste0("WT-1-", 1:3), paste0("Ao129-0-", 1:3), paste0("Ao129-1-", 1:3)),
            "Matrix time code 1 is mapped to 24 h using the GEO design."),
  make_spec("GSE215079", "GSE215079", "TPM", cfg_file("GSE215079"),
            c(paste0("AO_0_", 1:3), paste0("AO_12_", 1:3), paste0("AO_24_", 1:3), paste0("AO_6_", 1:3),
              paste0("age2_0_", 1:3), paste0("age2_12_", 1:3), paste0("age2_24_", 1:3), paste0("age2_6_", 1:3)),
            c(paste0("AO_0_", 1:3), paste0("AO_12_", 1:3), paste0("AO_24_", 1:3), paste0("AO_6_", 1:3),
              paste0("Q405_0_", 1:3), paste0("Q405_12_", 1:3), paste0("Q405_24_", 1:3), paste0("Q405_6_", 1:3)),
            "All samples are C. elegans induced; Q405 is retained as the mutant genotype label."),
  make_spec("GSE245014_Rab7", "GSE245014", "TPM", cfg_file("GSE245014_Rab7"),
            c(paste0("R7A_12_hypha", 1:3), paste0("R7A_36_hypha", 1:3), paste0("R7N_0_hypha", 1:3),
              paste0("WTA_12_hypha", 1:3), paste0("WTA_36_hypha", 1:3), paste0("WTN_0_hypha", 1:3)),
            c(paste0("Rab7_12_", 1:3), paste0("Rab7_36_", 1:3), paste0("Rab7_0_", 1:3),
              paste0("WT_12_", 1:3), paste0("WT_36_", 1:3), paste0("WT_0_", 1:3)),
            "Rab7 matrix uses R7A/R7N and WTA/WTN aliases; internal gene IDs are bridged from annotated tables."),
  make_spec("GSE245014_Vps35", "GSE245014", "TPM", cfg_file("GSE245014_Vps35"),
            c(paste0("AO_0h_", 1:3), paste0("AO_12h_", 1:3), paste0("AO_36h_", 1:3),
              paste0("V_0h_", 1:3), paste0("V_12h_", 1:3), paste0("V_36h_", 1:3)),
            c(paste0("WT_0_", 1:3), paste0("WT_12_", 1:3), paste0("WT_36_", 1:3),
              paste0("Vps35_0_", 1:3), paste0("Vps35_12_", 1:3), paste0("Vps35_36_", 1:3)),
            "Pseudo-XLS is a tab-delimited file; AO is WT and V is delta-Aovps35."),
  make_spec("GSE253759", "GSE253759", "TPM", cfg_file("GSE253759"),
            c(paste0("AO_0_", 1:3), paste0("AO_12_", 1:3), paste0("AO_24_", 1:3), paste0("AO_6_", 1:3),
              paste0("age2_0_", 1:3), paste0("age2_12_", 1:3), paste0("age2_24_", 1:3), paste0("age2_6_", 1:3)),
            c(paste0("AO_0_", 1:3), paste0("AO_12_", 1:3), paste0("AO_24_", 1:3), paste0("AO_6_", 1:3),
              paste0("Q405_0_", 1:3), paste0("Q405_12_", 1:3), paste0("Q405_24_", 1:3), paste0("Q405_6_", 1:3)),
            "Parallel-looking all-DEG table retained as a separate GEO record; do not pool as independent replicates."),
  make_spec("GSE301423_TPM", "GSE301423", "TPM", cfg_file("GSE301423"),
            c("WT_0h_1_tpm", "WT_0h_2_tpm", "WT_0h_3_tpm", "WT_24h_1_tpm", "WT_24h_2_tpm", "WT_24h_3_tpm",
              "ham4_0h_1_tpm", "ham4_0h_2_tpm", "ham4_0h_3_tpm", "ham4-24h-3_tpm", "ham4-24hh-2_tpm", "ham4_24h_1_tpm",
              "ham8-0h-1_tpm", "ham8-0h-2_tpm", "ham8-0h-3_tpm", "ham8-24h-1_tpm", "ham8-24h-2_tpm", "ham8-24h-3_tpm"),
            c("WT_0h_1_1", "WT_0h_2_1", "WT_0h_3_1", "WT_24h_1_1", "WT_24h_2_1", "WT_24h_3_1",
              "ham4_0h_1_1", "ham4_0h_2_1", "ham4_0h_3_1", "ham4-24h-3_1", "ham4-24h-2_1", "ham4_24h_1_1",
              "ham8-0h-1_1", "ham8-0h-2_1", "ham8-0h-3_1", "ham8-24h-1_1", "ham8-24h-2_1", "ham8-24h-3_1"),
            "Retains the two column typos ham4-24h-3 and ham4-24hh-2 with explicit GEO title mapping."),
  make_spec("GSE301423_counts", "GSE301423", "count", cfg_file("GSE301423"),
            c("WT_0h_1_count", "WT_0h_2_count", "WT_0h_3_count", "WT_24h_1_count", "WT_24h_2_count", "WT_24h_3_count",
              "ham4_0h_1_count", "ham4_0h_2_count", "ham4_0h_3_count", "ham4-24h-3_count", "ham4-24hh-2_count", "ham4_24h_1_count",
              "ham8-0h-1_count", "ham8-0h-2_count", "ham8-0h-3_count", "ham8-24h-1_count", "ham8-24h-2_count", "ham8-24h-3_count"),
            c("WT_0h_1_1", "WT_0h_2_1", "WT_0h_3_1", "WT_24h_1_1", "WT_24h_2_1", "WT_24h_3_1",
              "ham4_0h_1_1", "ham4_0h_2_1", "ham4_0h_3_1", "ham4-24h-3_1", "ham4-24h-2_1", "ham4_24h_1_1",
              "ham8-0h-1_1", "ham8-0h-2_1", "ham8-0h-3_1", "ham8-24h-1_1", "ham8-24h-2_1", "ham8-24h-3_1"),
            "Raw counts paired with the TPM columns above."),
  make_spec("GSE301452_FPKM", "GSE301452", "FPKM", cfg_file("GSE301452"),
            c("WT-0h-1_fpkm", "WT-0h-2_fpkm", "WT-0h-3_fpkm", "WT-24h-1_fpkm", "WT-24h-2_fpkm", "WT-24h-3_fpkm",
              "chk1-0h-1_fpkm", "chk1-0h-2_fpkm", "chk1-0h-3_fpkm", "chk1-24h-1_fpkm", "chk1-24h-2_fpkm", "chk1-24h-3_fpkm"),
            c("WT-0h-1_1", "WT-0h-2_1", "WT-0h-3_1", "WT-24h-1_1", "WT-24h-2_1", "WT-24h-3_1",
              "chk1-0h-1_1", "chk1-0h-2_1", "chk1-0h-3_1", "chk1-24h-1_1", "chk1-24h-2_1", "chk1-24h-3_1"),
            "FPKM values; the same file contains raw count columns."),
  make_spec("GSE301452_counts", "GSE301452", "count", cfg_file("GSE301452"),
            c("WT-0h-1_count", "WT-0h-2_count", "WT-0h-3_count", "WT-24h-1_count", "WT-24h-2_count", "WT-24h-3_count",
              "chk1-0h-1_count", "chk1-0h-2_count", "chk1-0h-3_count", "chk1-24h-1_count", "chk1-24h-2_count", "chk1-24h-3_count"),
            c("WT-0h-1_1", "WT-0h-2_1", "WT-0h-3_1", "WT-24h-1_1", "WT-24h-2_1", "WT-24h-3_1",
              "chk1-0h-1_1", "chk1-0h-2_1", "chk1-0h-3_1", "chk1-24h-1_1", "chk1-24h-2_1", "chk1-24h-3_1"),
            "Raw counts paired with the FPKM columns above."),
  make_spec("GSE303416_FPKM", "GSE303416", "FPKM", cfg_file("GSE303416"),
            c(paste0("A1-", 1:3, "_fpkm"), paste0("A2-", 1:3, "_fpkm"), paste0("STK1-", 1:3, "_fpkm"), paste0("STK2-", 1:3, "_fpkm")),
            c(paste0("Hal4b_0_", 1:3), paste0("Hal4b_12_", 1:3), paste0("WT_0_", 1:3), paste0("WT_12_", 1:3)),
            "A1/A2/STK1/STK2 are mapped in the GEO sample order (Hal4b 0/12 then WT 0/12); mapping confidence is recorded."),
  make_spec("GSE303416_counts", "GSE303416", "count", cfg_file("GSE303416"),
            c(paste0("A1-", 1:3, "_count"), paste0("A2-", 1:3, "_count"), paste0("STK1-", 1:3, "_count"), paste0("STK2-", 1:3, "_count")),
            c(paste0("Hal4b_0_", 1:3), paste0("Hal4b_12_", 1:3), paste0("WT_0_", 1:3), paste0("WT_12_", 1:3)),
            "Raw counts paired with the FPKM columns above."),
  make_spec("GSE306156_TPM", "GSE306156", "TPM", cfg_file("GSE306156"),
            c(paste0("AO-0h-", 1:3, "_tpm"), paste0("AO-18h-", 1:3, "_tpm"), paste0("Opy2-0h-", 1:3, "_tpm"), paste0("Opy2-18h-", 1:3, "_tpm")),
            c(paste0("AO_0h_", 1:3), paste0("AO_18h_", 1:3), paste0("Opy2_0h_", 1:3), paste0("Opy2_18h_", 1:3)),
            "WT/AoOpy2 mutant, 0/18 h nematode induction."),
  make_spec("GSE306156_counts", "GSE306156", "count", cfg_file("GSE306156"),
            c(paste0("AO-0h-", 1:3, "_count"), paste0("AO-18h-", 1:3, "_count"), paste0("Opy2-0h-", 1:3, "_count"), paste0("Opy2-18h-", 1:3, "_count")),
            c(paste0("AO_0h_", 1:3), paste0("AO_18h_", 1:3), paste0("Opy2_0h_", 1:3), paste0("Opy2_18h_", 1:3)),
            "Raw counts paired with the TPM columns above."),
  make_spec("GSE322548_TPM", "GSE322548", "TPM", cfg_file("GSE322548_TPM"),
            c(paste0("WT_", rep(c("0h", "12h", "24h"), each = 3), "_", rep(1:3, 3)), paste0("Hog1_", rep(c("0h", "12h", "24h"), each = 3), "_", rep(1:3, 3))),
            c(paste0("WT_", rep(c("0h", "12h", "24h"), each = 3), "_", rep(1:3, 3), "_1"), paste0("hog1_", rep(c("0h", "12h", "24h"), each = 3), "_", rep(1:3, 3), "_1")),
            "WT/Aohog1 KD, 0/12/24 h nematode induction."),
  make_spec("GSE322548_counts", "GSE322548", "count", cfg_file("GSE322548_counts"),
            c(paste0("AO_", rep(c("0", "12", "24"), each = 3), "_", rep(1:3, 3)), paste0("Hog1_", rep(c("0", "12", "24"), each = 3), "_", rep(1:3, 3))),
            c(paste0("WT_", rep(c("0h", "12h", "24h"), each = 3), "_", rep(1:3, 3), "_1"), paste0("hog1_", rep(c("0h", "12h", "24h"), each = 3), "_", rep(1:3, 3), "_1")),
            "Raw counts paired with the TPM columns above; AO is WT in the count file."),
  make_spec("GSE322740_FPKM", "GSE322740", "FPKM", cfg_file("GSE322740"),
            c(paste0("AO-", 1:3), paste0("M-", 1:3), paste0("D-", 1:3)),
            c(paste0("Fungal hyphae cells, AO-", 1:3), paste0("Fungal hyphae cells, M-", 1:3), paste0("Fungal hyphae cells, D-", 1:3)),
            "Day 5 aerial hyphae differentiation; no time-series contrast is implied.")
)

if (any(vapply(specs, function(s) is.na(s$file) || !file.exists(s$file), logical(1)))) stop("One or more configured files are missing")

normalize_title <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))
catalog_samples$sample_title_key <- normalize_title(catalog_samples$sample_title)

design_for <- function(gse, matrix_id, title, catalog_row) {
  t <- tolower(title)
  genotype <- "unknown"
  if (gse == "GSE210886") {
    genotype <- if (matrix_id == "GSE210886_PKA") if (grepl("^c_", t)) "delta_AoPkaC1" else "WT" else if (grepl("^acy", t)) "delta_AoAcy" else "WT"
  } else if (gse == "GSE212735") genotype <- if (grepl("^a11", t)) "delta_Aoatg11" else "WT"
  else if (gse == "GSE213542") genotype <- if (grepl("ao129", t)) "delta_AobZip129" else "WT"
  else if (gse %in% c("GSE215079", "GSE253759")) genotype <- if (grepl("q405|age2", t)) "delta_Aoage2" else "WT"
  else if (gse == "GSE245014") genotype <- if (grepl("rab7", t)) "delta_AoRab7" else if (grepl("vps35", t)) "delta_AoVps35" else "WT"
  else if (gse == "GSE301423") genotype <- if (grepl("ham4", t)) "delta_Aoham4" else if (grepl("ham8", t)) "delta_Aoham8" else "WT"
  else if (gse == "GSE301452") genotype <- if (grepl("chk1", t)) "delta_Aochk1" else "WT"
  else if (gse == "GSE303416") genotype <- if (grepl("hal4b", t)) "delta_Aohal4b" else "WT"
  else if (gse == "GSE306156") genotype <- if (grepl("opy2", t)) "delta_AoOpy2" else "WT"
  else if (gse == "GSE322548") genotype <- if (grepl("hog1", t)) "delta_Aohog1" else "WT"
  else if (gse == "GSE322740") genotype <- if (grepl("fungal hyphae cells, m-", t)) "delta_AoMdv1" else if (grepl("fungal hyphae cells, d-", t)) "delta_AoVps1" else "WT"
  time_h <- parse_time_hours(title)
  if (gse == "GSE213542" && grepl("-1-", title, fixed = TRUE)) time_h <- 24
  if (gse == "GSE212735" && grepl("_N_", title, fixed = TRUE)) time_h <- 12
  if (gse == "GSE322740") time_h <- NA_real_
  time_label <- if (gse == "GSE322740") "Day5" else if (is.finite(time_h)) paste0(time_h, "h") else "unknown"
  treatment <- if (gse == "GSE212735" && is.finite(time_h) && time_h > 0) "nitrogen_deficiency" else if (gse == "GSE322740") "aerial_hyphae_differentiation" else if (is.finite(time_h) && time_h > 0) "nematode_induction" else "baseline"
  if (gse %in% c("GSE215079", "GSE253759")) treatment <- "C_elegans_induction"
  if (gse == "GSE210886") treatment <- if (is.finite(time_h) && time_h > 0) "nematode_extract" else "baseline"
  list(genotype = genotype, time_h = time_h, time_label = time_label, treatment = treatment,
       baseline = isTRUE(!is.na(time_h) && time_h == 0), mapping_confidence = if (gse == "GSE303416") "order_based" else "title_based")
}

## Read each distinct file only once, before building a gene-ID bridge.
unique_files <- unique(vapply(specs, function(s) s$file, character(1)))
raw_tables <- setNames(lapply(unique_files, function(p) {
  tbl <- table_with_header(read_geo_table(p))
  attr(tbl, "source_file") <- p
  tbl
}), unique_files)
bridge <- make_gene_bridge(lapply(raw_tables, function(tbl) list(table = tbl)), candidates)
write_tsv(bridge, file.path(PROJ, "results", "expression", "gene_id_bridge.tsv"))

all_long <- list()
matrix_records <- list()
matrix_qc <- list()
sample_manifest <- list()
mapping_summary <- data.frame()

for (sp in specs) {
  cat("parsing", sp$matrix_id, "\n")
  tbl <- raw_tables[[sp$file]]
  missing_cols <- setdiff(sp$cols, colnames(tbl))
  if (length(missing_cols)) stop("Missing columns in ", sp$matrix_id, ": ", paste(missing_cols, collapse = ", "))
  rowmap <- map_rows_to_candidates(tbl, candidates, bridge)
  idcols <- grep("gene.?id|^id$", colnames(tbl), ignore.case = TRUE, value = TRUE)
  namecols <- grep("gene.?name|symbol", colnames(tbl), ignore.case = TRUE, value = TRUE)
  if (!length(idcols)) idcols <- colnames(tbl)[1]
  raw_id <- as.character(tbl[[idcols[1]]])
  raw_name <- if (length(namecols)) as.character(tbl[[namecols[1]]]) else raw_id
  rowmap$raw_gene_id <- raw_id
  rowmap$raw_gene_name <- raw_name
  rowmap$row_index <- seq_len(nrow(tbl))
  mapped_counts <- table(rowmap$candidate_id)
  duplicate_candidates <- names(mapped_counts)[mapped_counts > 1 & nzchar(names(mapped_counts))]
  values <- sapply(sp$cols, function(nm) as_numeric_safely(tbl[[nm]]))
  values <- as.matrix(values)
  colnames(values) <- sp$cols
  rownames(values) <- make.unique(raw_id)
  if (sp$unit == "count") values <- round(pmax(values, 0))
  matrix_records[[sp$matrix_id]] <- list(matrix_id = sp$matrix_id, gse = sp$gse, unit = sp$unit,
                                          file = sp$file, columns = sp$cols, sample_titles = sp$titles,
                                          values = values, rowmap = rowmap, notes = sp$notes)
  mapped <- rowmap$candidate_id
  candidate_row <- split(which(!is.na(mapped)), mapped[!is.na(mapped)])
  spec_sample_rows <- list()
  for (j in seq_along(sp$cols)) {
    title <- sp$titles[j]
    key <- normalize_title(title)
    hit <- which(catalog_samples$gse == sp$gse & catalog_samples$sample_title_key == key)
    if (!length(hit)) stop("GEO sample title not found: ", sp$gse, " / ", title)
    if (length(hit) > 1) stop("Ambiguous GEO sample title: ", sp$gse, " / ", title)
    cr <- catalog_samples[hit, ]
    ds <- design_for(sp$gse, sp$matrix_id, title, cr)
    repnum <- suppressWarnings(as.integer(sub(".*[^0-9]([0-9]+)(?:_[0-9]+)?$", "\\1", gsub("_(?:tpm|fpkm|count)$", "", sp$cols[j], ignore.case = TRUE))))
    if (!is.finite(repnum)) repnum <- j
    sample_key <- paste(sp$gse, sp$matrix_id, cr$gsm, sp$unit, sep = "::")
    sm <- data.frame(gse = sp$gse, matrix_id = sp$matrix_id, gsm = cr$gsm, sample_title = title,
                     matrix_column = sp$cols[j], sample_key = sample_key, unit = sp$unit,
                     genotype = ds$genotype, time_h = ds$time_h, time_label = ds$time_label,
                     treatment = ds$treatment, baseline = ds$baseline, replicate = repnum,
                     mapping_confidence = ds$mapping_confidence, source_file = sp$file,
                     source_sha256 = file_sha256(sp$file), stringsAsFactors = FALSE)
    spec_sample_rows[[j]] <- sm
    rows <- candidate_row
    for (cid in candidates$sequence_id) {
      ri <- rows[[cid]]
      if (is.null(ri)) {
        val <- NA_real_; rid <- NA_character_; rnm <- NA_character_; mm <- "not_in_expression_matrix"
      } else {
        val <- if (length(ri) == 1) values[ri, j] else mean(values[ri, j], na.rm = TRUE)
        if (!is.finite(val)) val <- NA_real_
        rid <- paste(unique(raw_id[ri]), collapse = ";")
        rnm <- paste(unique(raw_name[ri]), collapse = ";")
        mm <- paste(unique(rowmap$mapping_method[ri]), collapse = ";")
      }
      all_long[[length(all_long) + 1L]] <- data.frame(candidate_id = cid,
          candidate_core_id = candidates$candidate_core_id[match(cid, candidates$sequence_id)],
          gse = sp$gse, matrix_id = sp$matrix_id, gsm = cr$gsm, sample_title = title,
          matrix_column = sp$cols[j], sample_key = sample_key, unit = sp$unit,
          genotype = ds$genotype, time_h = ds$time_h, time_label = ds$time_label,
          treatment = ds$treatment, baseline = ds$baseline, replicate = repnum,
          value = val, raw_gene_id = rid, raw_gene_name = rnm, mapping_method = mm,
          detected_in_matrix = !is.na(val), expressed_gt_zero = is.finite(val) && val > 0,
          source_file = sp$file, stringsAsFactors = FALSE)
    }
  }
  sample_manifest[[sp$matrix_id]] <- do.call(rbind, spec_sample_rows)
  mapped_unique <- unique(na.omit(mapped))
  mapping_summary <- rbind(mapping_summary, data.frame(matrix_id = sp$matrix_id, gse = sp$gse,
      unit = sp$unit, rows = nrow(tbl), columns = ncol(tbl), expression_columns = length(sp$cols),
      candidates_mapped = length(mapped_unique), candidates_missing = nrow(candidates) - length(mapped_unique),
      duplicate_candidate_rows = length(duplicate_candidates), direct_rows = sum(rowmap$mapping_method == "direct_id_column"),
      bridge_rows = sum(rowmap$mapping_method == "gene_id_bridge"), nonnegative_fraction = mean(values >= 0, na.rm = TRUE),
      numeric_fraction = mean(is.finite(values)), file_sha256 = file_sha256(sp$file), notes = sp$notes,
      stringsAsFactors = FALSE))
}

long <- rbindlist(all_long, fill = TRUE)
samples <- rbindlist(sample_manifest, fill = TRUE)
setorder(long, candidate_id, gse, matrix_id, unit, gsm)
write_tsv(long, file.path(PROJ, "results", "expression", "candidate_expression_long.tsv"))
long_norm <- long[unit != "count"]
write_tsv(long_norm, file.path(PROJ, "results", "expression", "candidate_expression_normalized_long.tsv"))
long_wide <- dcast(long, candidate_id + candidate_core_id ~ sample_key, value.var = "value")
write_tsv(long_wide, file.path(PROJ, "results", "expression", "candidate_expression_wide.tsv"))
write_tsv(samples, file.path(PROJ, "results", "expression", "candidate_sample_metadata.tsv"))
write_tsv(mapping_summary, file.path(PROJ, "results", "expression", "matrix_parse_qc.tsv"))

candidate_map <- merge(candidates, unique(long[, .(candidate_id, matrix_id, mapping_method, detected_in_matrix)], by = c("candidate_id", "matrix_id")),
                       by.x = "sequence_id", by.y = "candidate_id", all.x = TRUE, allow.cartesian = TRUE)
write_tsv(candidate_map, file.path(PROJ, "results", "expression", "candidate_id_mapping.tsv"))
saveRDS(matrix_records, file.path(PROJ, "results", "expression", "parsed_matrices.rds"), compress = "xz")
write_tsv(data.frame(path = c(file.path(PROJ, "results", "expression", "candidate_expression_long.tsv"),
                              file.path(PROJ, "results", "expression", "candidate_expression_wide.tsv"),
                              file.path(PROJ, "results", "expression", "parsed_matrices.rds")),
                       sha256 = vapply(c(file.path(PROJ, "results", "expression", "candidate_expression_long.tsv"),
                                         file.path(PROJ, "results", "expression", "candidate_expression_wide.tsv"),
                                         file.path(PROJ, "results", "expression", "parsed_matrices.rds")), file_sha256, character(1))),
          file.path(PROJ, "checksums", "step02_outputs_sha256.tsv"))

cat("spec_count:", length(specs), "\n")
cat("sample_metric_rows:", nrow(samples), "\n")
cat("candidate_expression_rows:", nrow(long), "\n")
cat("candidate_count:", length(unique(long$candidate_id)), "\n")
cat("missing_EVM_rows:", sum(long$mapping_method == "not_in_expression_matrix"), "\n")
print(mapping_summary[, c("matrix_id", "unit", "expression_columns", "candidates_mapped", "candidates_missing", "direct_rows", "bridge_rows")])
finish_log()
