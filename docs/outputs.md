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
                                     # (both in the curate step)
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
- `best_hit_kingdom` — Bacteria/Archaea/Virus/Eukaryote/Unknown (from the FoldSeek
  lineage, else the taxid, else the organism name; `all_hits.csv` has it per hit as
  `kingdom` + `kingdom_source`)

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
| `foldseek` | FoldSeek-only: named from the FoldSeek hits, no informative Phold call |
| `phold` | Phold-only: Phold call, no FoldSeek support |
| `pharokka` | non-hypothetical, kept from Pharokka |
| `manual_review` / `curation_suggestion` | needs_review rows (reviewer-edited / rule engine's suggestion accepted as-is) |
| `structural-uncharacterised` | a CONFIDENT/GOOD fold match to an **unnamed** target (written by `07_integrate`). `final_product` keeps the target's header (often a DUF), so the fold is recorded — but the FUNCTION is unknown. **Do not count these as annotated.** |
| `no_informative_hit` | still hypothetical: hits, but none naming a function |
| `no_hit` | still hypothetical: no hit |

## `04_output/updated_prophages.gb`
GenBank with curated `/product`, `/function`, and `/note` provenance per CDS.
Built by one of three paths (update Pharokka GBK / concat per-prophage GBKs /
build from coords for the pre-CDS sub-mode). Not produced in pangenome sub-mode
(no genome coordinates).

## `04_output/pangenome_reference.{csv,xlsx,faa}` (pangenome sub-mode)
Built after curation in pangenome mode, from the id map
(`<run>/split/pangenome/pangenome_id_map.csv`). Rejoins the final table (keyed on the
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
  `phynteny_input.fasta` (the FASTA it wrote).

## `04_output/curation/review_suggested.csv`
Genes flagged for review, with the proposed name (`final_annotation`) and the reason
(`curation_explanation`). To change a call, add a row to `input/overrides.tsv` and
re-run (only curation re-runs; the review flag is cleared).

## Interactive report (`04_output/report.html`)

One self-contained HTML file (data embedded, opens offline, can be emailed): genome
map per prophage coloured by function (genes in input order when the run has no
GenBank), a searchable/filterable table, and per gene the evidence (Phold, FoldSeek
best hit + top-3 with links to AlphaFold DB / RCSB, curation reason, phynteny).
"Export visible rows" saves the current selection as CSV. Built by the Snakemake
`report` rule, or for any run dir: `python scripts/09_report.py --run_dir <run_dir>`.
