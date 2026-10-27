# Database setup (local FoldSeek search mode)

WebAPI search mode (`02w`) needs no local FoldSeek databases; Pharokka and Phold
still need theirs (and the optional taxdump).

## What you need

All installed by `phagefactor databases`, controlled by flags so you only
download what you need (`--all` = everything except the afdb50 index):

| DB | Purpose | Approx size | Flag |
|---|---|---|---|
| **pharokka_db** | Pharokka annotation (genome mode) | ~1–5 GB | `--pharokka` |
| **phold_db** | PHold structural search (otherwise phold's own default location) | ~2 GB | `--phold` |
| **pdb100** | experimental PDB structures | ~4 GB | `--pdb100` |
| **afdb-swissprot** | curated SwissProt AFDB predictions | ~1.5 GB | `--swissprot` |
| **afdb50** | bacterial/archaeal/viral AFDB depth (UniProt50-minimal) | ~40–60 GB | `--afdb50` |
| **afdb50 index** | precomputed `.idx` + `_ss.idx`, loaded into RAM by the search (256 GB node) | ~350 GB | `--index` |
| **taxdump** | NCBI taxonomy: hit domain when FoldSeek gives no lineage (webapi mode, older runs). Optional | ~250 MB | `--taxdump` |

```bash
phagefactor databases /data/phagefactor_db --pharokka --phold --pdb100 --swissprot --afdb50 --taxdump
sbatch --mem=200G --time=04:00:00 --wrap "phagefactor databases --index"   # afdb50 index, fat node
phagefactor check                                                          # what is installed
```

The directory is recorded in `~/.config/phagefactor/settings.yaml` and used by every
project (`<dir>/pharokka_db`, `<dir>/phold_db`, `<dir>/foldseek_dbs`); later calls
without a directory add to it. A project can override it in `config/config.yaml →
databases` (or with `PHAGEFACTOR_DB_ROOT`, `FOLDSEEK_DB_ROOT`, ...); confirm the
per-DB relative paths in `databases.foldseek_local_dbs`.

> **Critical:** each value must be the FoldSeek **DB root**, NOT the `*_ca` file.
> FoldSeek auto-discovers `_ss` / `_h` from the root name. The `_ca` (Cα
> coordinate) files are not used and can be deleted to save space.

## Taxonomy & the taxon filter

Each FoldSeek hit gets a domain (`kingdom`) and its source (`kingdom_source`), from
the first available of (`scripts/lib/taxonomy.py`):

1. **lineage**: FoldSeek's `taxlineage` (local search; afdb50, afdb-swissprot and
   pdb100 have taxonomy).
2. **taxid**: walked up the NCBI taxdump (web-API mode, FFAL tables, older runs).
   `phagefactor databases --taxdump`; optional.
3. **name**: organism-name keywords (`lexicon.py`), which miss many eukaryotes.

Otherwise the domain is Unknown, never treated as eukaryotic. The 02 log counts
each source.

`databases.foldseek_taxon_filter` is a per-DB keep-list of domain taxids (2 =
Bacteria, 2157 = Archaea, 10239 = Viruses, 2759 = Eukaryota); afdb50 keeps
`2,2157,10239`, dropping hits of a known other domain at parse time.

## RAM

The afdb50 search needs a large-memory node. Two things keep it there and not higher:

1. **Precomputed index** (recommended): build once with
   `foldseek createindex <db> <tmp> --index-exclude 2 --threads 16`.
   This produces both `<db>.idx` (AA) and `<db>_ss.idx` (3Di) — FoldSeek then
   searches in low-RAM mode instead of building the k-mer table in memory.
   (`phagefactor databases --index` does this.)
2. **`--sort-by-structure-bits 0`**: when supported, FoldSeek does not load Cα into
   RAM (~151 GB → ~35 GB for afdb50). `02d` auto-detects support and adds it.

> **Don't try to shrink the RAM request below the index size.** The default
> `foldseek_search` step requests 256 GB (`resources:` in `workflow.yaml`) because the afdb50 3Di index is
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
