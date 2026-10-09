# Modes

Three independent choices:

| axis | options | set by |
|---|---|---|
| **Input mode** | genome · protein · pangenome | auto-detected, or `INPUT_MODE=` |
| **Runner** | HPC (SLURM) · workstation | `submit_all.sh` · `run_local.sh` |
| **Search backend** | on-disk DBs · web server · none | `SEARCH_MODE=` |

They combine freely: a workstation run can use on-disk DBs, and an HPC run is
always on-disk (compute nodes have no internet).

## Input mode (auto-detected by `detect_input_type.py`)

Detection on `input/fasta/`: any sequence > 5 kb → **genome**; otherwise if
sequences are amino-acid (≥2 % protein-only residues) → **protein**; short
ATCGN-only contigs default to genome. Result cached in `input/.input_type`.
Override with `INPUT_MODE=genome|protein` or `--force`.

### Genome mode (mode 1)
One nucleotide FASTA per prophage in `input/fasta/<NAME>.fasta`, names in
`input/prophage_list.txt`.
```
00c pharokka_array → 01 phold_array → 01c merge + 01d merge_3di → 02 → 03 → 04 → 05 → 06 → 07
```
- Pharokka uses **Phanotate** (catches small phage ORFs; do not pass `--fast`).
- Phold runs GPU+autotune (`phold_use_gpu: 1`) or falls back to CPU.
- The FoldSeek **target set = original Pharokka hypotheticals** (not Phold's
  post-rescue product) — this preserves the merged-evidence rescues even as newer
  Phold versions annotate more proteins up front.

### Protein mode (mode 2 / bulk CDS)
Bulk protein FASTA; **all sequences treated as hypothetical candidates**.
```
00p split (50/batch) → 01p phold-proteins array → 01p merge → 02 → 03 → 04 → 05 → 06 → 07
```
- Batch size is fixed (`protein_batch_size: 50`) for reproducibility; last batch
  may be partial. CPU by default (`phold_proteins_use_gpu: 0`).
- **Pre-CDS-computed sub-mode**: if you supply a source-genome FASTA + prophage
  window coordinates (+ a rich gene-metadata CSV with absolute coords),
  `05_build_output.py` (part of `steps/04_curate.sh`) rebuilds a GenBank so
  genome-style outputs work even though there was no Pharokka run (used for
  the citro phiNP/phiSM eggNOG inputs).

### Pangenome sub-mode (bulk CDS, no genome coordinates)

A clustered pan-genome / pan-proteome (e.g. a Roary representative-sequence set:
one gene per record, many short in-frame CDS in a single FASTA) runs as **protein
mode**, with a translate step in front and a rejoin step at the end:

```
00n translate (nt CDS -> protein, clean ids)  →  [normal protein chain: 00p → 01p → 02 → 03 → 04]  →  08 rejoin (group + seq + annotation)
```

**Just drop the nucleotide FASTA in `input/fasta/` and run** (leave `INPUT_MODE`
unset, or set `INPUT_MODE=pangenome`). Auto-detection returns `pangenome`;
`submit_all.sh` then translates it and continues as protein mode automatically.

- **`scripts/00n_translate_pangenome.py`** — protein mode is *protein-in* (it does
  not translate). Translates the nucleotide CDS (default codon table 11), renames
  every record to a clean `pang_NNNNN` id (dotted/pipe ids like `x.gff.1` /
  `fig|..` can be mangled by phold/foldseek id handling), and writes
  `input/pangenome_id_map.csv` (`clean_id,orig_id,group,nt_len,aa_len`) plus
  `input/.input_type = protein`. `submit_all.sh` runs it automatically for
  `pangenome` mode and moves the original nucleotide file to `input/pangenome_nt/`
  so the protein split does not re-process it.
- **`detect_input_type.py`** returns the `pangenome` verdict for a bulk nucleotide
  CDS set (≥ `PANGENOME_MIN_RECORDS`, default 20, short nucleotide-only records),
  instead of silently choosing genome mode and Phanotate-recalling every record.
- **`00p_split_protein_batches.py`** guards protein mode against nucleotide input:
  if the FASTA is DNA it stops before any job with a pointer to the pangenome path.
- **No `prophage_windows.csv` / `genome/`** → step 05 skips the GenBank and there
  is **no phynteny / synteny step** (they need per-gene genome coordinates). The
  `final_annotations_table.{csv,xlsx}` is still produced in full.
- **`scripts/08_pangenome_reference.py`** — submitted automatically after curation
  when an id map is present (`steps/08_pangenome_reference.sh`). Rejoins the final
  table (keyed on `pang_NNNNN`) to the id map + protein FASTA and emits the
  **homology-reference set**: `pangenome_reference.{csv,xlsx}` (one row per group:
  `group_id, representative_locus, aa_sequence, final_product, short_name,
  annotation_source, agreement`) and `pangenome_reference.faa` (header
  `>group_id product`, ready for `diamond makedb` / `mmseqs createdb`). Search a
  larger protein collection against it to transfer annotations by homology; each
  hit's `sseqid` is the `group_id` to join back to the CSV.

## Runner

| | command | notes |
|---|---|---|
| **HPC (SLURM)** | `bash submit_all.sh` | one job per step, arrays for pharokka/phold |
| **Workstation** | `bash run_local.sh` | same steps, sequential, no queue |

`run_local.sh` is for testing and small jobs — full-depth annotation needs
afdb50, which needs cluster storage and memory.

## Search backend (step 02)

### On-disk databases (`02d_foldseek_3di.py`) — `SEARCH_MODE=local`, default
Searches the Phold 3Di tokens against installed DBs (PDB100, AFDB-SwissProt,
AFDB50). Builds the query DB via the official ProstT5 `tsv2db` recipe
(no `_ca` file). Per-DB taxonomy + taxon filtering, per-residue masking. Needs the
DBs on disk — see [`databases.md`](databases.md).

### Web server (`02w_foldseek_webapi.py`) — `SEARCH_MODE=webapi`
Submits the AA FASTA to `search.foldseek.com`; ProstT5 + search run server-side.
Ticket/poll with resubmit-on-timeout. Output: `webapi_hits.m8`, bridged via
`webm8_to_best_hit()` through the same `foldseek_scoring.build_best_and_top3`
path as the on-disk backend to produce `best_hit.csv`/`top3.csv` — step 03 runs
identically regardless of which search mode produced them. **Experimental.** The server searches structures, not sequences, so queries are
folded with ESMFold first (`scripts/lib/esmfold_cache.py`, cached by sequence).
The public ESMFold service is intermittent and caps at ~400 aa, which excludes
large structural proteins.

On a workstation use `SEARCH_MODE=local` with pdb100 + afdb-swissprot (~5.5 GB);
`python scripts/tools/check_databases.py` lists what is installed.

## Quick selector

| You have… | Input | Runner | Backend |
|---|---|---|---|
| Per-prophage genomes, cluster + DBs | genome | `submit_all.sh` | `local` |
| Per-prophage genomes, workstation, small DBs | genome | `run_local.sh` | `local` |
| Bulk CDS / pangenome proteins, cluster | protein | `submit_all.sh` | `local` |
| eggNOG-annotated source genome (citro) | protein (pre-CDS sub-mode) | either | `local` |
| No databases at all | any | `run_local.sh` | `webapi` |
