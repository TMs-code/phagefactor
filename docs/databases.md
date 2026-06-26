# Database setup (local FoldSeek search mode)

WebAPI search mode (`02w`) needs **no** local databases — skip this page if you
only use WebAPI.

## What you need

| DB | Purpose | Approx size | Install |
|---|---|---|---|
| **pharokka_db** | Pharokka annotation (genome mode) | ~1–5 GB | `install_databases.py -o <path>` |
| **pdb100** | experimental PDB structures | ~4 GB | `foldseek databases pdb100 <path>/pdb100_db /tmp` |
| **afdb-swissprot** | curated SwissProt AFDB predictions | ~1.5 GB | `foldseek databases afdb-swissprot <path>/afdb_swissprot_db /tmp` |
| **afdb50** (or BakTFold AFDBClusters) | bacterial AFDB depth | ~40–60 GB (afdb50) | `foldseek databases afdb50 <path>/afdb50_db /tmp` |

Point `config/config.yaml → databases.foldseek_db_root` (or export
`FOLDSEEK_DB_ROOT`) at the parent folder, and confirm the per-DB relative paths in
`databases.foldseek_local_dbs`.

> **Critical:** each value must be the FoldSeek **DB root**, NOT the `*_ca` file.
> FoldSeek auto-discovers `_ss` / `_ca` / `_h` from the root name. Pointing at
> `<root>_ca` silently breaks 3Di alignment.

## Taxonomy & the taxon filter

DBs with a `<root>_taxonomy` companion file get per-DB taxon filtering
(`databases.foldseek_taxon_filter`, e.g. afdb50 → bacteria+archaea+viruses only).
DBs without taxonomy (pdb100, afdb-swissprot here) keep all hits; the eukaryotic
filters in step 03/scoring handle them by description. `00c_add_taxonomy_baktfold.py`
attaches taxonomy to a BakTFold DB if needed.

## RAM and the no-index / low-RAM option

Large DBs (afdb50) can OOM during search. Two levers:

1. **Precomputed index** (recommended): build once with
   `foldseek createindex <db> <tmp> --index-exclude 2 --threads 16`.
   This produces both `<db>.idx` (AA) and `<db>_ss.idx` (3Di) — FoldSeek then
   searches in low-RAM mode instead of building the k-mer table in memory.
   (See `00e_createindex_afdb50.sh`.)
2. **`--sort-by-structure-bits 0`**: when supported, FoldSeek does not load Cα into
   RAM (~151 GB → ~35 GB for afdb50). `02d` auto-detects support and adds it.

`02d` logs a full DB-file audit before each search (which `.idx` files exist, sizes,
incomplete-index detection) so you can diagnose OOM/index issues from the log.

## Note on the query DB

`02d` builds the query DB with the official ProstT5 `foldseek tsv2db` recipe (AA +
3Di + headers, no `_ca`). This works on FoldSeek 10.x without a ProstT5-enabled
build, as long as no Cα-derived output column is requested (the format string is
sequence-only by design).
