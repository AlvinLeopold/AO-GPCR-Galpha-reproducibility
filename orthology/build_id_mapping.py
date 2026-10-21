#!/usr/bin/env python3
"""Build protein/gene ID mapping tables from an OrthoFinder result directory."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


RESULTS = Path(".")
OUT = Path(__file__).resolve().parent
SPECIES = {
    "0": ("Arthrobotrys_oligospora_AOL", "Arthrobotrys_oligospora_AOL.fasta"),
    "1": ("Arthrobotrys_oligospora_TWF", "Arthrobotrys_oligospora_TWF.fasta"),
}
AOL = SPECIES["0"][0]
TWF = SPECIES["1"][0]


def split_ids(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def read_sequence_ids() -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    path = RESULTS / "WorkingDirectory" / "SequenceIDs.txt"
    for line in path.read_text().splitlines():
        if not line.strip() or ":" not in line:
            continue
        internal_id, header = line.split(":", 1)
        internal_id, header = internal_id.strip(), header.strip()
        species_idx = internal_id.split("_", 1)[0]
        species, fasta = SPECIES[species_idx]
        protein_id = header.split()[0]
        match = re.search(r"(?:^|\s)gene=([^\s]+)", header)
        gene_id = match.group(1) if match else protein_id
        records[protein_id] = {
            "internal_id": internal_id,
            "species": species,
            "fasta_file": fasta,
            "protein_id": protein_id,
            "gene_id": gene_id,
            "original_header": header,
        }
    return records


def read_orthogroups(path: Path) -> dict[str, dict[str, list[str]]]:
    groups: dict[str, dict[str, list[str]]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            og = row["Orthogroup"]
            groups[og] = {species: split_ids(row.get(species, "")) for species in (AOL, TWF)}
    return groups


def read_inferred_orthologs() -> dict[str, dict[str, list[str]]]:
    """Read AOL-to-TWF calls and return members/partners by orthogroup."""
    path = RESULTS / "Orthologues" / f"{AOL}.tsv"
    calls: dict[str, dict[str, list[str]]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            og = row["Orthogroup"]
            calls[og] = {
                AOL: split_ids(row[AOL]),
                TWF: split_ids(row["Orthologs"]),
            }
    return calls


def relation(aol_n: int, twf_n: int) -> str:
    if aol_n == 1 and twf_n == 1:
        return "one-to-one"
    if aol_n == 1 and twf_n > 1:
        return "one-to-many"
    if aol_n > 1 and twf_n == 1:
        return "many-to-one"
    if aol_n > 1 and twf_n > 1:
        return "many-to-many"
    return "species-specific"


def main() -> None:
    global RESULTS, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True,
                        help="Path to the OrthoFinder results directory")
    parser.add_argument("--out", type=Path, default=Path.cwd(),
                        help="Output directory for four TSV mapping tables")
    args = parser.parse_args()
    RESULTS = args.results.resolve()
    OUT = args.out.resolve()
    required = [
        RESULTS / "WorkingDirectory/SequenceIDs.txt",
        RESULTS / "Orthogroups/Orthogroups.tsv",
        RESULTS / "Orthogroups/Orthogroups_UnassignedGenes.tsv",
        RESULTS / "Orthologues/Arthrobotrys_oligospora_AOL.tsv",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        parser.error("Missing OrthoFinder inputs: " + ", ".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)
    sequence = read_sequence_ids()
    assigned = read_orthogroups(RESULTS / "Orthogroups" / "Orthogroups.tsv")
    unassigned = read_orthogroups(RESULTS / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv")
    inferred = read_inferred_orthologs()

    groups: dict[str, dict[str, list[str]]] = {}
    status: dict[str, str] = {}
    protein_to_og: dict[str, str] = {}
    for og, members in assigned.items():
        groups[og] = members
        for protein_id in members[AOL] + members[TWF]:
            status[protein_id] = "assigned"
            protein_to_og[protein_id] = og
    for og, members in unassigned.items():
        groups[og] = members
        for protein_id in members[AOL] + members[TWF]:
            status[protein_id] = "unassigned"
            protein_to_og[protein_id] = og

    # Inferred partners indexed by protein ID, preserving all one-to-many links.
    partners: dict[str, list[str]] = defaultdict(list)
    for members in inferred.values():
        for aol_id in members[AOL]:
            partners[aol_id].extend(members[TWF])
        for twf_id in members[TWF]:
            partners[twf_id].extend(members[AOL])
    partners = {k: sorted(set(v)) for k, v in partners.items()}

    # One row per input protein, including proteins that OrthoFinder left unassigned.
    protein_rows: list[dict[str, str | int]] = []
    for protein_id, record in sequence.items():
        og = protein_to_og.get(protein_id, "")
        members = groups.get(og, {AOL: [], TWF: []})
        species = record["species"]
        other_species = TWF if species == AOL else AOL
        all_other = members[other_species]
        inferred_other = partners.get(protein_id, [])
        if status.get(protein_id) == "unassigned":
            call_status = "unassigned"
        elif inferred_other:
            call_status = "inferred_ortholog"
        elif all_other:
            call_status = "shared_orthogroup_only"
        else:
            call_status = "species_specific"
        protein_rows.append(
            {
                **record,
                "orthogroup": og,
                "assignment_status": status.get(protein_id, "not_found"),
                "orthogroup_relation": relation(len(members[AOL]), len(members[TWF])),
                "orthogroup_other_species_proteins": ", ".join(all_other),
                "inferred_ortholog_partners": ", ".join(inferred_other),
                "inferred_partner_count": len(inferred_other),
                "ortholog_call_status": call_status,
            }
        )

    protein_rows.sort(key=lambda x: (x["species"], x["internal_id"]))
    protein_header = [
        "internal_id", "species", "fasta_file", "protein_id", "gene_id", "original_header",
        "orthogroup", "assignment_status", "orthogroup_relation", "orthogroup_other_species_proteins",
        "inferred_ortholog_partners", "inferred_partner_count", "ortholog_call_status",
    ]
    with (OUT / "protein_id_mapping.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=protein_header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(protein_rows)

    # Explicit pairwise rows from OrthoFinder's Orthologues output.
    pair_rows: list[dict[str, str]] = []
    for og, members in inferred.items():
        aol_ids, twf_ids = members[AOL], members[TWF]
        rel = relation(len(aol_ids), len(twf_ids))
        for aol_id in aol_ids:
            for twf_id in twf_ids:
                pair_rows.append(
                    {
                        "orthogroup": og,
                        "relation_type": rel,
                        "AOL_protein_id": aol_id,
                        "AOL_gene_id": sequence.get(aol_id, {}).get("gene_id", ""),
                        "TWF_protein_id": twf_id,
                        "TWF_gene_id": sequence.get(twf_id, {}).get("gene_id", ""),
                    }
                )
    pair_rows.sort(key=lambda x: (x["orthogroup"], x["AOL_protein_id"], x["TWF_protein_id"]))
    pair_header = ["orthogroup", "relation_type", "AOL_protein_id", "AOL_gene_id", "TWF_protein_id", "TWF_gene_id"]
    with (OUT / "AOL_TWF_ortholog_pairs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_header, delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows)

    # One row per orthogroup, covering all assigned groups (including groups
    # for which OrthoFinder did not emit an explicit Orthologues call).
    group_rows: list[dict[str, str | int]] = []
    pair_counts = Counter(row["orthogroup"] for row in pair_rows)
    inferred_ogs = set(inferred)
    for og, members in assigned.items():
        aol_ids, twf_ids = members[AOL], members[TWF]
        group_rows.append(
            {
                "orthogroup": og,
                "AOL_count": len(aol_ids),
                "TWF_count": len(twf_ids),
                "relation_type": relation(len(aol_ids), len(twf_ids)),
                "AOL_protein_ids": ", ".join(aol_ids),
                "TWF_protein_ids": ", ".join(twf_ids),
                "explicit_orthologues_call": "yes" if og in inferred_ogs else "no",
                "explicit_ortholog_pair_count": pair_counts[og],
            }
        )
    group_rows.sort(key=lambda x: x["orthogroup"])
    group_header = [
        "orthogroup", "AOL_count", "TWF_count", "relation_type", "AOL_protein_ids", "TWF_protein_ids",
        "explicit_orthologues_call", "explicit_ortholog_pair_count",
    ]
    with (OUT / "AOL_TWF_orthogroup_mapping.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=group_header, delimiter="\t")
        writer.writeheader()
        writer.writerows(group_rows)

    # Compact summary for reproducibility and quick QC.
    summary = [
        ("results_directory", str(RESULTS)),
        ("AOL_proteins", sum(r["species"] == AOL for r in protein_rows)),
        ("TWF_proteins", sum(r["species"] == TWF for r in protein_rows)),
        ("total_proteins", len(protein_rows)),
        ("assigned_proteins", sum(r["assignment_status"] == "assigned" for r in protein_rows)),
        ("unassigned_proteins", sum(r["assignment_status"] == "unassigned" for r in protein_rows)),
        ("assigned_orthogroups", len(assigned)),
        ("shared_orthogroups", sum(bool(v[AOL]) and bool(v[TWF]) for v in assigned.values())),
        ("inferred_ortholog_orthogroups", len(inferred)),
        ("inferred_ortholog_pairs", len(pair_rows)),
        ("assigned_orthogroups_with_explicit_call", len(inferred_ogs)),
    ]
    counts = Counter(row["relation_type"] for row in pair_rows)
    for key in ("one-to-one", "one-to-many", "many-to-one", "many-to-many"):
        summary.append((f"inferred_pairs_{key}", counts[key]))
    with (OUT / "mapping_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)

    print(f"Wrote {len(protein_rows)} protein rows and {len(pair_rows)} inferred ortholog pairs to {OUT}")


if __name__ == "__main__":
    main()
