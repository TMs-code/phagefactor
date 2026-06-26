# Example 2 — protein / pre-annotated-CDS mode (2 Citrobacter prophages)

Use this mode when the CDS are **already defined** (e.g. a host genome annotation
that was quantified by RNA-seq), so re-calling genes would break the read mapping.
Pharokka is skipped; you supply the proteins plus their coordinates/annotations, and
phold + custom FoldSeek + the FACT curation layer annotate them. A GenBank is rebuilt
from the source-genome slice so phynteny + C1/Cro synteny still run.

Worked example: the two spontaneously-inducible *Citrobacter rodentium* prophages
**phiNP** and **phiSM** (pre-annotated via eggNOG from the host genome).

## Input (shipped, ready to run)
```
protein_2citro/input/
├── .input_type                      # one line: protein
├── fasta/
│   ├── phiNP_proteins.faa           # 105 proteins; header id = locus_tag (gene_id_fa)
│   ├── phiSM_proteins.faa           # 68 proteins
│   └── prophage_windows.csv         # prophage,genome_start,genome_stop
├── gene_metadata_rich.csv           # the join table, 173 rows (schema in ../README.md)
└── genome/
    └── source_genome.fasta          # Citrobacter rodentium host genome — for the GBK rebuild
```

## Run it
```bash
cp -r example_data/protein_2citro/input/* input/
bash submit_all.sh                   # auto-detects "protein"
```
`gene_metadata_rich.csv` is the join-key table: `gene_id_fa` matches each protein FASTA
header; `start/end/strand/contig` give the host-genome coordinates used to rebuild the
per-prophage GenBank (so phynteny + C1/Cro synteny still run). Full column list in `../README.md`.

