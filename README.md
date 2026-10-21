# AO GPCR–Gα reproducibility code

Reusable scripts and notebooks for **“Genome-Wide GPCR Mapping Links AoGPR200 and AoGPRPR8 to Trap Formation and Distinct Gα Responses in *Arthrobotrys oligospora*”** (IJMS manuscript ijms-4611842).

## Data

The associated [Zenodo dataset](https://doi.org/10.5281/zenodo.23524073) contains three archives:

1. `AO_GPCR_transcriptomic_source_and_processed_2026-10-21.tar.gz` — a complete runnable package with 15 unchanged GEO supplementary expression/count source files from 12 series, sample metadata, candidate and processed ranking tables, Figure 2 source data, Supplementary Data S1/S2, compact structural scores/inputs, scripts, provenance manifests, and integrity checker. These GEO files are **processed matrices, not raw FASTQ reads**.
2. `AO_GPCR_Galpha_AF2_complete_2026-10-21.tar.zst` — the complete original LocalColabFold/AlphaFold2 output directory, including 267 complex jobs, 801 ranked model structures/scores, PAE, MSAs, run scripts and logs.
3. `AOL_TWF_OrthoFinder_ID_Mapping_complete_2026-10-21.tar.gz` — both input proteomes, complete OrthoFinder output and original mapping tables/script.

The Zenodo record also includes SHA-256 checksums and a data README. Original GEO accessions and file-level SHA-256 values appear in `R2_computation/inputs/atcc_gpcr_analysis/metadata/source_file_manifest.tsv` within the first archive.

## Reproduce

For a complete rerun, download the first Zenodo archive, then extract it into a separate working directory:

```bash
mkdir AO_GPCR_reproducibility
tar -xzf AO_GPCR_transcriptomic_source_and_processed_2026-10-21.tar.gz
cd AO_GPCR_reproducibility
python3 verify_package.py
./reproduce_r2.sh
```

The upstream 12-series R pipeline starts at `R2_computation/inputs/atcc_gpcr_analysis/run_all.sh`; its dependencies are listed in `environment/README.md` inside the data archive. The bundled original run logs document the historical execution. The packaging machine validated the R2 Python analyses and Figure 2 regeneration; it did **not** rerun the full upstream R/Bioconductor pipeline or LocalColabFold.

To regenerate the AOL–TWF mapping tables from the third archive:

```bash
python3 orthology/build_id_mapping.py \
  --results /path/to/02_orthofinder_complete_results \
  --out /path/to/reproduced_mapping
```

The three regenerated mapping tables matched the originals byte for byte in a local test. The summary table also matched after excluding its machine-specific `results_directory` row.

## License and scope

Repository code is licensed under MIT. The Zenodo dataset is licensed under CC BY 4.0. Cite the original GEO accessions for the source expression/count matrices. Experimental BRET, trap-count and gel/microscopy originals are outside this computational release.
