# Output reference

## Folder layout

```
04_output/
├── final_annotations_table.csv     # deliverable table
├── final_annotations_table.xlsx    # colour-coded
├── updated_prophages.gb            # curated GenBank (not produced in pangenome sub-mode)
├── pangenome_reference.csv         # pangenome sub-mode only: one row per group
├── pangenome_reference.xlsx        #   (see below)
├── pangenome_reference.faa         #   DIAMOND/MMseqs-ready reference FASTA
└── curation/
    ├── curated_annotations.csv     # written by 04_curate_annotations.py
    └── review_suggested.csv        # written by 05_build_output.py
                                     # (both run in one job, steps/04_curate.sh)
05_phynteny/
├── phynteny.tsv  phynteny.gbk  phynteny_predictions.csv
├── final_annotations_integrated.csv  # phynteny merged into final_function
├── final_with_synteny.gb             # GBK with synteny/phynteny notes
└── run/
    ├── phynteny_run.log            # phynteny_transformer's verbose log
    └── phynteny_input.fasta        # the FASTA phynteny wrote (relocated here)
```

## `04_output/final_annotations_table.{csv,xlsx}`

One row per CDS. XLSX is colour-coded (resolved = green, still-hypothetical =
orange). Key columns:

**Identity & answer**
- `prophage`, `locus_tag`, `aa_length`
- `short_name` — concise label (see naming rules)
- `final_product`, `final_function` — the curated call
- `was_hypothetical` — was this a Pharokka hypothetical (i.e. a pipeline target)?
- `annotation_source` — evidence category (see below)
- `best_hit_kingdom` — Bacteria/Archaea/Virus/Eukaryote/Unknown

**Pharokka**: `pharokka_product`, `pharokka_function`
**Phold**: `phold_product`, `phold_function_cat`, `phold_confidence`, `phold_phrog`,
`accession_phrog`, `phold_evalue`
**Curation**: `agreement` (from step 03, carried verbatim), `curation_action` (which
rule resolved the gene), `curation_suggestion` / `curation_explanation` (the rule
engine's proposed call and its stated reasoning — always deterministic, rule-based;
no LLM/AI is involved in generating these)
**Custom FoldSeek**: `foldseek_description`, `foldseek_accession`, `foldseek_db`,
`foldseek_taxname`, `foldseek_confidence`, `foldseek_score`, `foldseek_evalue`,
`foldseek_pident`, `foldseek_qcov_frac`, `foldseek_partial_match`,
`foldseek_same_host`, `foldseek_top3` (a.k.a. `fs_top3`),
`fs_top3_kingdoms` — the kingdom (Bacteria/Archaea/Viruses/Eukaryota/Unknown) of
each of the FoldSeek top-3 hits, positionally parallel to `fs_top3`
**Provenance**: `note` (`pipeline=phagefactor; source=…; …; original_pharokka=hypothetical protein`)

### Evidence categories (`annotation_source`)
| Value | Meaning |
|---|---|
| `both agree` | Phold and custom-FoldSeek concur (strong/partial) |
| `merged` | combined into one more-specific call |
| `foldseek` → **structural-only** | FoldSeek structural homology, no Phold corroboration |
| `phold` → **phold-only** | Phold call, no FoldSeek support |
| `pharokka` | non-hypothetical, kept from Pharokka |
| `manual_review` / `curation_suggestion` | needs_review rows (reviewer-edited / rule engine's suggestion accepted as-is) |
| `structural-uncharacterised` | a CONFIDENT/GOOD fold match to an **unnamed** target (written by `07_integrate`). `final_product` keeps the target's header (often a DUF), so the fold is recorded — but the FUNCTION is unknown. **Do not count these as annotated.** |
| `no_hit` | still hypothetical |

> Reader-facing rename: surface `foldseek` as **structural-only** in figures/tables.

## `04_output/updated_prophages.gb`
GenBank with curated `/product`, `/function`, and `/note` provenance per CDS.
Built by one of three paths (update Pharokka GBK / concat per-prophage GBKs /
build from coords for the pre-CDS sub-mode). Not produced in pangenome sub-mode
(no genome coordinates).

## `04_output/pangenome_reference.{csv,xlsx,faa}` (pangenome sub-mode)
Built by `steps/08_pangenome_reference.sh` after curation, when an id map
(`input/pangenome_id_map.csv`) is present. Rejoins the final table (keyed on the
clean `pang_NNNNN` ids) back to the biological group names + sequences.

- `.csv` / `.xlsx` — one row per pan-genome group: `group_id`,
  `representative_locus`, `clean_id`, `aa_length`, `aa_sequence`, `final_product`,
  `short_name`, `annotation_source`, `agreement` (plus any evidence columns present).
- `.faa` — the same representatives as a FASTA, header `>group_id  final_product`.
  Ready for `diamond makedb` / `mmseqs createdb`; a homology hit's `sseqid` **is**
  the `group_id`, so annotation transfer to a larger protein set is a plain join
  back to the CSV.

## Phynteny step (`05_phynteny/`)
- `05_phynteny/phynteny_predictions.csv` — synteny category + probability (≥0.8).
- `05_phynteny/final_annotations_integrated.csv` — phynteny category merged into
  `final_function` (phynteny priority on confident gaps) + `synteny_hint` (C1/Cro,
  lysis cassette, integration module).
- `05_phynteny/final_with_synteny.gb` — GBK with synteny/phynteny notes.
- `05_phynteny/run/` — `phynteny_run.log` (phynteny_transformer's verbose log) and
  `phynteny_input.fasta` (the FASTA it wrote), relocated/renamed there by
  `steps/05_phynteny.sh`.

## `04_output/curation/review_suggested.csv`
Genes flagged for manual review (written by `05_build_output.py`, part of
`steps/04_curate.sh`). Read `curation_explanation` for the rule engine's stated
reasoning, accept/edit `final_annotation`, then re-run. Pre-filled with
`curation_suggestion`. (`04_curate_annotations.py`, the other half of the same
job, writes the companion `04_output/curation/curated_annotations.csv`.)
