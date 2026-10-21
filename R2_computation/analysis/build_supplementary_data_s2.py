#!/usr/bin/env python3
"""Build a compact, reproducible Supplementary Data S2 archive."""

from pathlib import Path
from hashlib import sha256
from shutil import copy2, rmtree
from zipfile import ZipFile, ZIP_DEFLATED
import csv

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUP = ROOT / "supplementary"
PACKAGE = SUP / "Supplementary_Data_S2_Computational_Sensitivity"
ZIP = SUP / "Supplementary_Data_S2_Computational_Sensitivity.zip"

if PACKAGE.exists():
    rmtree(PACKAGE)
PACKAGE.mkdir(parents=True)

rank = pd.read_csv(ROOT / "analysis/transcriptome_rank_sensitivity.tsv", sep="\t")
rank.to_csv(SUP / "Supplementary_Table_S5A_Rank_Sensitivity.tsv", sep="\t", index=False)
support = pd.read_csv(ROOT / "analysis/transcriptome_exclusion_sensitivity.tsv", sep="\t")
support.to_csv(SUP / "Supplementary_Table_S5B_Exclusion_Sensitivity.tsv", sep="\t", index=False)

files = [
    "inputs/selection_criterion.tsv",
    "inputs/selection_regulator_effects.tsv",
    "inputs/Figure2_SourceData_PanelC.tsv",
    "inputs/tree_annotation_metadata.tsv",
    "inputs/af2_all_model_scores.tsv",
    "inputs/af2_focal_seed101_scores.tsv",
    "inputs/focal_seed101_run_command.txt",
    "inputs/submitted_structure_summary.tsv",
    "inputs/submitted_coupling_selectivity.tsv",
    "inputs/config.json",
    "inputs/run_stage1.sh",
    "inputs/score_complexes.py",
    "inputs/atcc_gpcr_analysis/results/effects/candidate_effects.tsv",
    "inputs/atcc_gpcr_analysis/results/integration/candidate_evidence_by_gse.tsv",
    "inputs/atcc_gpcr_analysis/results/figure_source_data/Fig2a_WT_induction_atlas.tsv",
    "inputs/atcc_gpcr_analysis/results/design/series_design_manifest.tsv",
    "analysis/reanalyse_computation.py",
    "analysis/analyse_focal_seed.py",
    "analysis/make_figure2_r2.R",
    "analysis/make_af2_sensitivity_figure.R",
    "analysis/transcriptome_regulator_effects.tsv",
    "analysis/transcriptome_rank_sensitivity.tsv",
    "analysis/transcriptome_exclusion_sensitivity.tsv",
    "analysis/transcriptome_leave_one_series_out.tsv",
    "analysis/low_expression_pseudocount_sensitivity.tsv",
    "analysis/af2_all_model_scores_clean.tsv",
    "analysis/af2_pair_all_model_summary.tsv",
    "analysis/af2_receptor_rank_robustness.tsv",
    "analysis/af2_focal_pair_sensitivity.tsv",
    "analysis/af2_focal_two_seed_scores.tsv",
    "analysis/af2_focal_seed_sensitivity.tsv",
    "analysis/af2_focal_seed_pair_ranges.tsv",
    "analysis/af2_focal_seed_summary.json",
    "analysis/summary.json",
    "analysis/figure2_9series/Figure2_SourceData_PanelA.tsv",
    "analysis/figure2_9series/Figure2_SourceData_PanelB.tsv",
    "analysis/figure2_9series/Figure2_SourceData_PanelC.tsv",
    "analysis/figure2_9series/Figure2_SourceData_PanelD.tsv",
    "analysis/figure2_9series/Figure2_SourceData_StudyDesign.tsv",
    "analysis/af2_figure/Supplementary_Figure_AF2_model_variability.pdf",
]
files += [str(p.relative_to(ROOT)) for p in sorted((ROOT / 'inputs/focal_a3m').glob('*.a3m'))]
assert len([f for f in files if f.endswith('.a3m')]) == 6
files += [str(p.relative_to(ROOT)) for p in sorted((ROOT / 'inputs/focal_structures').rglob('*.pdb'))]
assert len([f for f in files if f.endswith('.pdb')]) == 12

for rel in files:
    src = ROOT / rel
    if not src.is_file():
        raise FileNotFoundError(src)
    dest = PACKAGE / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    copy2(src, dest)

for name in ("Supplementary_Table_S5A_Rank_Sensitivity.tsv", "Supplementary_Table_S5B_Exclusion_Sensitivity.tsv", "Supplementary_Table_S6_AF2_Seed_Sensitivity.tsv"):
    copy2(SUP / name, PACKAGE / name)

readme = """Supplementary Data S2  Computational sensitivity analyses

Contents
- inputs/: archived processed contrast tables, Figure 2 inputs, all 801 original AF2 model ipTM scores, 18 focal seed-101 scores, six focal A3M inputs, 12 top-ranked focal PDB structures, original run configuration, and scoring script.
- analysis/: executable reanalysis and figure scripts, full source tables, the nine-series Figure 2 panel data, Supplementary Figure S3 PDF, and focal seed comparison.
- Supplementary_Table_S5A_Rank_Sensitivity.tsv and S5B_Exclusion_Sensitivity.tsv: the two parts of Table S5.
- Supplementary_Table_S6_AF2_Seed_Sensitivity.tsv: focal seed-0 versus seed-101 scores and rankings.

Reproduce the tabular analyses from this archive root with:
  python3 analysis/reanalyse_computation.py
  python3 analysis/analyse_focal_seed.py

Requires pandas and NumPy. The code reads only files inside this archive and writes to analysis/.
The R figure scripts use base R. To regenerate the revised nine-series Figure 2, pass the seven paths documented at the top of analysis/make_figure2_r2.R.

Method provenance
- AF2 screen: LocalColabFold 1.5.5, commit e8ebd9a612b9da5fe12bbbc09513d5c52a9af633; models 1-3, one seed (0), no templates, 3 recycles.
- Focal six-pair sensitivity rerun: same six archived A3M files, models 1-3 and 3 recycles, random seed 101. The exact command is in inputs/focal_seed101_run_command.txt.
- All-model score extraction was read from 801 archived rank-score JSONs on the authors' workstation. The focal top-ranked PDBs from both seeds are included; the complete 801-pair structure collection and full PAE JSONs are not contained in this compact archive.
- GSE253759 and GSE303416 are excluded from the revised Figure 2 WT atlas. The original regulator-effect ranking uses GSE210886, GSE306156 and GSE322548, so exclusion affects WT-response support rather than those pathway-effect values.

Limitations
- Structural model weights and the two stochastic seeds are not independent experimental replicates or a functional coupling benchmark.
- The scoring thresholds 0.02 and 0.05 are uncalibrated and are not used as biological selectivity classes.
- BRET and trap raw measurements are outside this archive; no claims about their independent biological replication are made here.
"""
(PACKAGE / "README.txt").write_text(readme)

manifest = []
for path in sorted(PACKAGE.rglob("*")):
    if path.is_file():
        manifest.append((path.relative_to(PACKAGE).as_posix(), sha256(path.read_bytes()).hexdigest(), path.stat().st_size))
with (PACKAGE / "SHA256SUMS.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["relative_path", "sha256", "bytes"])
    w.writerows(manifest)

with ZipFile(ZIP, "w", compression=ZIP_DEFLATED, compresslevel=6) as z:
    for path in sorted(PACKAGE.rglob("*")):
        if path.is_file():
            z.write(path, path.relative_to(PACKAGE))
print(ZIP)
print("files", len(manifest) + 1, "zip_bytes", ZIP.stat().st_size)
