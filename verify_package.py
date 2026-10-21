#!/usr/bin/env python3
"""Validate reviewer-facing source files and archive inventory using stdlib only."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "R2_computation/inputs"
UPSTREAM = INPUTS / "atcc_gpcr_analysis"
MANIFEST = ROOT / "checksums/package_sha256.tsv"
QA = ROOT / "QA_report.md"
EXCLUDED = {MANIFEST, QA}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    errors: list[str] = []
    original = UPSTREAM / "metadata/source_file_manifest.tsv"
    with original.open(newline="") as handle:
        geo_rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in geo_rows:
        path = (UPSTREAM / row["file"]).resolve()
        if not path.is_relative_to(INPUTS.resolve()) or not path.is_file():
            errors.append(f"GEO source missing or outside package: {row['source_id']}")
            continue
        if path.stat().st_size != int(row["file_size_bytes"]):
            errors.append(f"GEO size mismatch: {row['source_id']}")
        if sha256(path) != row["sha256"]:
            errors.append(f"GEO SHA-256 mismatch: {row['source_id']}")

    if len(geo_rows) != 15 or len({row["gse"] for row in geo_rows}) != 12:
        errors.append("Expected 15 GEO files from 12 GSE accessions")

    with (UPSTREAM / "checksums/input_sha256.tsv").open(newline="") as handle:
        input_rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in input_rows:
        path = INPUTS / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            errors.append(f"Missing or size-mismatched upstream input: {row['input_role']}")
        elif sha256(path) != row["sha256"]:
            errors.append(f"SHA-256 mismatch for upstream input: {row['input_role']}")

    scores = INPUTS / "af2_all_model_scores.tsv"
    with scores.open(newline="") as handle:
        score_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(score_rows) != 801 or len({row["job"] for row in score_rows}) != 267:
        errors.append("Expected 801 AF2 model scores in 267 jobs")

    s1 = ROOT / "supplementary/S1_R2"
    s1_old = ROOT / "supplementary/Supplementary_Data_S1_qRTPCR_Figure2.zip"
    with ZipFile(s1_old) as archive:
        if archive.testzip() is not None:
            errors.append("Prior S1 archive failed CRC check")
        for name in archive.namelist():
            if name.startswith("qRT_PCR/") and (s1 / name).read_bytes() != archive.read(name):
                errors.append(f"R2 S1 qRT-PCR source changed: {name}")
    for panel in (s1 / "Figure2_processed").glob("*.tsv"):
        original_panel = ROOT / "R2_computation/analysis/figure2_9series" / panel.name
        if not original_panel.is_file() or sha256(panel) != sha256(original_panel):
            errors.append(f"R2 Figure 2 source mismatch: {panel.name}")

    for archive_path in sorted((ROOT / "supplementary").glob("*.zip")):
        with ZipFile(archive_path) as archive:
            bad = archive.testzip()
            if bad is not None:
                errors.append(f"ZIP CRC failure: {archive_path.name}: {bad}")
    with ZipFile(ROOT / "supplementary/Supplementary_Data_S1_R2_candidate.zip") as archive:
        expected_names = {p.relative_to(s1).as_posix() for p in s1.rglob("*") if p.is_file()}
        if set(archive.namelist()) != expected_names:
            errors.append("R2 S1 ZIP inventory differs from the browseable S1_R2 directory")
        for name in expected_names:
            if archive.read(name) != (s1 / name).read_bytes():
                errors.append(f"R2 S1 ZIP content differs: {name}")

    files = sorted(
        p for p in ROOT.rglob("*")
        if p.is_file() and p not in EXCLUDED
        and "reproduced" not in p.relative_to(ROOT).parts
        and ".git" not in p.relative_to(ROOT).parts
    )
    oversized = [p.relative_to(ROOT).as_posix() for p in files if p.stat().st_size >= 100_000_000]
    if oversized:
        errors.append(f"Files exceed GitHub's 100 MB file limit: {oversized}")

    current = {p.relative_to(ROOT).as_posix(): (p.stat().st_size, sha256(p)) for p in files}
    if MANIFEST.exists() and "--refresh" not in sys.argv:
        with MANIFEST.open(newline="") as handle:
            expected = {
                row["path"]: (int(row["size_bytes"]), row["sha256"])
                for row in csv.DictReader(handle, delimiter="\t")
            }
        if expected != current:
            errors.append("Package manifest differs; inspect edits, then run verify_package.py --refresh")
    else:
        MANIFEST.parent.mkdir(exist_ok=True)
        with MANIFEST.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["path", "size_bytes", "sha256"])
            for name, (size, digest) in current.items():
                writer.writerow([name, size, digest])

    status = "PASS" if not errors else "FAIL"
    QA.write_text(
        "# Package QA\n\n"
        f"Status: **{status}**\n\n"
        f"Files checked: **{len(files)}**; GEO source files: **{len(geo_rows)}** "
        f"across **{len({row['gse'] for row in geo_rows})}** GSEs; "
        f"AF2 model score rows: **{len(score_rows)}**.\n\n"
        "R2 sensitivity Python scripts and Figure 2 regeneration were validated separately; "
        "the full upstream R/Bioconductor pipeline and 11 GB AF2 stage were not rerun on this Mac.\n\n"
        + ("Errors:\n\n" + "\n".join(f"- {e}" for e in errors) + "\n" if errors else "No integrity errors detected.\n")
    )
    print(f"{status}: {len(files)} files; {len(geo_rows)} GEO files; {len(score_rows)} AF2 rows")
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
