# Database setup (local FoldSeek search mode)

WebAPI search mode (`02w`) needs **no** local databases — skip this page if you
only use WebAPI.

## What you need

All four installed by one script, `steps/00b_setup_databases.sh`, controlled by
flags so you only download what you need (or `--all` for the first three):

| DB | Purpose | Approx size | Flag |
|---|---|---|---|
| **pharokka_db** | Pharokka annotation (genome mode) | ~1–5 GB | `--pharokka` |
| **pdb100** | experimental PDB structures | ~4 GB | `--pdb100` |
| **afdb-swissprot** | curated SwissProt AFDB predictions | ~1.5 GB | `--swissprot` |
| **afdb50** | bacterial/archaeal/viral AFDB depth (UniProt50-minimal) | ~40–60 GB | `--afdb50` |

```bash
bash steps/00b_setup_databases.sh --pharokka --pdb100 --swissprot --afdb50
```

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
filters in step 03/scoring handle them by description instead.

## RAM and the no-index / low-RAM option

Large DBs (afdb50) can OOM during search. Two levers:

1. **Precomputed index** (recommended): build once with
   `foldseek createindex <db> <tmp> --index-exclude 2 --threads 16`.
   This produces both `<db>.idx` (AA) and `<db>_ss.idx` (3Di) — FoldSeek then
   searches in low-RAM mode instead of building the k-mer table in memory.
   (See `steps/00b_setup_databases.sh --index`.)
2. **`--sort-by-structure-bits 0`**: when supported, FoldSeek does not load Cα into
   RAM (~151 GB → ~35 GB for afdb50). `02d` auto-detects support and adds it.

`02d` logs a full DB-file audit before each search (which `.idx` files exist, sizes,
incomplete-index detection) so you can diagnose OOM/index issues from the log.

> **Don't try to shrink the RAM request below the index size.** The default
> `steps/02_foldseek_3di.sh` requests 256 GB because the afdb50 3Di index is
> ~196 GB and is loaded fully into RAM for the prefilter. Switching to mmap
> (`--db-load-mode 2`) to schedule on a smaller node does **not** help here —
> at less-than-index-size RAM the index thrashes from disk and a search that
> takes minutes can take well over an hour. If you don't have a node with
> enough RAM, use WebAPI search mode instead of trying to shrink this.

## Note on the query DB

`02d` builds the query DB with the official ProstT5 `foldseek tsv2db` recipe (AA +
3Di + headers, no `_ca`). This works on FoldSeek 10.x without a ProstT5-enabled
build, as long as no Cα-derived output column is requested (the format string is
sequence-only by design).
