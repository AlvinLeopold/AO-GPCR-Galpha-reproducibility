# Software environment

The original upstream workstation logs record **R 4.4.1 on Ubuntu 24.10**. The package versions recorded in `R2_computation/inputs/atcc_gpcr_analysis/metadata/software_versions.tsv` include readxl 1.5.0, data.table 1.18.4, Biostrings 2.74.1, DESeq2 1.46.0, ggplot2 4.0.3, pheatmap 1.0.13 and jsonlite 2.0.0. Full `sessionInfo()` output for each upstream step is retained under `logs/`.

The upstream scripts import **readxl, data.table, Biostrings, digest, limma, DESeq2, ggplot2, patchwork, ragg, svglite, pheatmap and scales**. WGCNA, GSVA and fgsea paths are conditional and their availability/status is written to the output status tables. Install matching CRAN/Bioconductor versions before running `atcc_gpcr_analysis/run_all.sh`.

The R2 Python scripts require pandas and NumPy. Packaging validation used **Python with pandas 2.2.3 and NumPy 2.3.5**. The revised Figure 2 renderer uses base R. AF2 prediction was run on the workstation with **LocalColabFold 1.5.5**, commit `e8ebd9a612b9da5fe12bbbc09513d5c52a9af633`; model/seed/recycle settings are recorded in the package's structural input files and S2 README.

The code is distributed as analysis scripts and archived results rather than as an `renv` lockfile. Reproducing the full upstream result on another operating system may require resolving compatible binary packages. The original run logs and byte-matched source matrices provide the reference for that check.
