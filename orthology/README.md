# AOL–TWF OrthoFinder ID mapping

The complete OrthoFinder results, both input proteomes, and the original mapping
tables are archived in the associated Zenodo dataset. The script here regenerates
the four mapping tables from the `02_orthofinder_complete_results` directory in
that archive without relying on a workstation-specific path.

```bash
python3 orthology/build_id_mapping.py \
  --results /path/to/02_orthofinder_complete_results \
  --out /path/to/reproduced_mapping
```

The script uses only the Python standard library. It reads OrthoFinder's
`WorkingDirectory/SequenceIDs.txt`, `Orthogroups/Orthogroups.tsv`,
`Orthogroups/Orthogroups_UnassignedGenes.tsv`, and
`Orthologues/Arthrobotrys_oligospora_AOL.tsv`. It preserves explicit
one-to-many and many-to-many ortholog calls, and does not equate membership in
the same orthogroup with an explicit ortholog call.

The supplied Zenodo archive is the original source supplied for this revision.
The script was adapted only to accept `--results` and `--out` paths and validate
the required inputs. The original unmodified script is retained inside the
Zenodo archive for provenance.
