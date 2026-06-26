# example_data — worked examples + required input structure per mode

phageFACTor runs in three input modes. The **input files and folder layout differ
per mode** — this is the most important thing for a new user, so each is spelled out
below with a small worked example. Auto‑detection (`detect_input_type.py`) reads
`input/fasta/` and picks genome vs protein; pangenome is explicit.

```
example_data/
├── genome_2campy/     # Example 1: genome mode  (2 Campylobacter prophages) — shipped, ready to run
├── protein_2citro/    # Example 2: protein / pre‑annotated‑CDS mode (2 Citrobacter prophages: phiNP, phiSM)
├── run_smoke_test.sh  # end‑to‑end check on Example 1 (WebAPI mode, no local DB)
└── pangenome_TODO/    # Example 3: pangenome mode (bulk CDS, no genome coords) — to add
```

---
## Example 1 — GENOME mode (the simplest; 2 Campylobacter prophages)
Input = just the nucleotide genome FASTAs. Pharokka (Phanotate) calls the CDS, phold
+ FoldSeek annotate. **Minimal input:**
```
input/
├── fasta/
│   ├── CCIE1.fasta            # one nucleotide FASTA per prophage (>name on the header)
│   └── CCIE2.fasta
└── prophage_list.txt          # the prophage basenames, one per line (CCIE1 / CCIE2)
```
`.input_type` is written automatically (= `genome`). Run: `bash submit_all.sh`.

---
## Example 2 — PROTEIN / pre‑annotated‑CDS mode (2 Citrobacter prophages: phiNP, phiSM)
Use when the CDS are **predefined** (e.g. host annotation quantified by RNAseq) — no
pharokka. You provide the proteins + their coordinates/annotations. (Input files are
not committed yet — see `protein_2citro/README.md`.) **Required input:**
```
input/
├── .input_type                       # contains the word: protein
├── fasta/
│   ├── <Prophage>_proteins.faa       # AA FASTA per prophage; header = gene id (locus_tag)
│   └── prophage_windows.csv          # prophage,genome_start,genome_stop,contig
├── gene_metadata_rich.csv            # the join table (schema below)
└── genome/
    └── source_genome.fasta           # host genome (multi‑contig OK) — for the GBK
```
`gene_metadata_rich.csv` columns (one row per CDS):
`prophage, contig, gene_id_gff, gene_id_fa, start, end, strand, aa_length,
gff_product, emapper_desc, COG_category, preferred_name, is_hypo_gff,
is_hypo_emapper, is_hypothetical`
- `gene_id_fa` must equal the protein FASTA header id (join key).
- `start/end` are absolute host‑genome coords (windows can span the origin); `contig`
  is optional for a single-chromosome host (as in the citro example).
- `is_hypothetical` (= hypo in the host annotation AND in eggNOG) defines the FoldSeek
  target set.
- You build this table once from your host annotation (GFF + proteins.fa + eggNOG +
  the prophage coordinate windows); `protein_2citro/` is a worked example of the result.
Run: `bash submit_all.sh` (auto‑detects protein).

---
## Example 3 — PANGENOME mode (TODO)
Bulk CDS with **no per‑gene genome coordinates** (e.g. a clustered pan‑proteome).
Same as protein mode minus `prophage_windows.csv` + `genome/` → **no GBK / no phynteny /
no synteny step**; FoldSeek + curation still produce the annotation table. To be added
once Tan has a pangenome dataset.

---
## Expected output (any mode)
`04_output/final_annotations_table.csv` (one row per CDS) + `04_output/curation/review_suggested.csv`
(the few genes flagged for a human look) + `05_phynteny/final_annotations_integrated.csv`
after the phynteny step. `run_smoke_test.sh` checks the table's schema and annotated
fraction; it runs WebAPI mode (no local DB) so a check needs no big database.

> Status: both `genome_2campy/` and `protein_2citro/` (phiNP + phiSM, 173 proteins) are
> shipped and runnable. After a run you may commit a small `final_annotations_table.csv`
> subset under an `expected/` folder as an illustrative reference.
