# Modes

Three independent choices:

| axis | options | set by |
|---|---|---|
| **Input mode** | genome · protein · pangenome | auto-detected, or `input_mode` |
| **Runner** | HPC (SLURM) · workstation | `phagefactor run --slurm` · `phagefactor run` |
| **Search backend** | on-disk DBs · web server · none | `search_mode` |

`input_mode` and `search_mode` are set in `config/workflow.yaml` or for one launch
with `--config search_mode=webapi`.

They combine freely. With `--slurm`, `search_mode: local` runs the FoldSeek search
as a compute-node job; `search_mode: webapi` runs it on the login node, since
compute nodes usually have no internet (the other steps still go to compute nodes).

## Input mode (auto-detected by `detect_input_type.py`)

Detection on `input/fasta/`: any sequence > 5 kb → **genome**; otherwise if
sequences are amino-acid (≥2 % protein-only residues) → **protein**; short
ATCGN-only contigs default to genome. Detection is redone at every launch and
writes nothing into the input dir. Override with `input_mode: genome|protein|pangenome`.

### Genome mode (mode 1)
One nucleotide FASTA per prophage in `input/fasta/<NAME>.fasta`, names in
`input/prophage_list.txt`.
```
00c pharokka_array → 01 phold_array → 01c merge + 01d merge_3di → 02 → 03 → 04 → 05 → 06 → 07 → 09 report
```
- Pharokka uses **Phanotate** (catches small phage ORFs; do not pass `--fast`). FASTA
  headers are reduced to the sequence id first (Phold rejects colons in the GenBank).
- Phold runs on CPU by default; `phold_gpu: true` for GPU + autotune.
- The FoldSeek **target set = original Pharokka hypotheticals** (not Phold's
  post-rescue product) — this preserves the merged-evidence rescues even as newer
  Phold versions annotate more proteins up front.

### Protein mode (mode 2 / bulk CDS)
Bulk protein FASTA; **all sequences treated as hypothetical candidates**.
```
00p split (50/batch) → 01p phold-proteins array → 01p merge → 02 → 03 → 04 → 05 → 06 → 07 → 09 report
```
- Batch size is fixed (`protein_batch_size: 50`) for reproducibility; last batch
  may be partial. CPU by default (`phold_proteins_gpu: true` for GPU).
- **Pre-CDS-computed sub-mode**: if you supply a source-genome FASTA + prophage
  window coordinates (+ a rich gene-metadata CSV with absolute coords),
  `05_build_output.py` (part of the curate step) rebuilds a GenBank so
  genome-style outputs work even though there was no Pharokka run (used for
  the citro phiNP/phiSM eggNOG inputs).

### Pangenome sub-mode (bulk CDS, no genome coordinates)

A clustered pan-genome / pan-proteome (e.g. a Roary representative-sequence set:
one gene per record, many short in-frame CDS in a single FASTA) runs as **protein
mode**, with a translate step in front and a rejoin step at the end:

```
00n translate (nt CDS -> protein, clean ids)  →  [normal protein chain: 00p → 01p → 02 → 03 → 04]  →  08 rejoin (group + seq + annotation)
```

**Just drop the nucleotide FASTA in `input/fasta/` and run** (leave `input_mode`
on `auto`, or set it to `pangenome`). Auto-detection returns `pangenome`; the
workflow translates it and continues as protein mode.

- **`scripts/00n_translate_pangenome.py`** — protein mode is *protein-in* (it does
  not translate). Translates the nucleotide CDS (default codon table 11), renames
  every record to a clean `pang_NNNNN` id (dotted/pipe ids like `x.gff.1` /
  `fig|..` can be mangled by phold/foldseek id handling), and writes
  `<run>/split/pangenome/pangenome_id_map.csv` (`clean_id,orig_id,group,nt_len,aa_len`).
  The translated FASTA also stays in the run dir; the input dir is not touched.
- **`detect_input_type.py`** returns the `pangenome` verdict for a bulk nucleotide
  CDS set (≥ `PANGENOME_MIN_RECORDS`, default 20, short nucleotide-only records),
  instead of silently choosing genome mode and Phanotate-recalling every record.
- **`00p_split_protein_batches.py`** guards protein mode against nucleotide input:
  if the FASTA is DNA it stops before any job with a pointer to the pangenome path.
- **No `prophage_windows.csv` / `genome/`** → step 05 skips the GenBank and there
  is **no phynteny / synteny step** (they need per-gene genome coordinates). The
  `final_annotations_table.{csv,xlsx}` is still produced in full.
- **`scripts/08_pangenome_reference.py`** — submitted automatically after curation
  in pangenome mode. Rejoins the final
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
| **HPC (SLURM)** | `phagefactor run --slurm` | one job per step, one per prophage / protein batch for pharokka and phold |
| **Workstation** | `phagefactor run` | same steps on this machine, in parallel up to `--cores` |

The workstation runner is for testing and small jobs: full-depth annotation needs
afdb50, which needs cluster storage and memory.

## Search backend (step 02)

### On-disk databases (`02d_foldseek_3di.py`) — `search_mode: local`, default
Searches the Phold 3Di tokens against installed DBs (PDB100, AFDB-SwissProt,
AFDB50). Builds the query DB via the official ProstT5 `tsv2db` recipe
(no `_ca` file). Per-hit lineage, per-DB taxon filtering, per-residue masking. Needs
the DBs on disk — see [`databases.md`](databases.md).

### Web server (`02w_foldseek_webapi.py`) — `search_mode: webapi`
Searches pdb100, afdb50 and afdb-swissprot on `search.foldseek.com`: no local FoldSeek
database or index. Pharokka and PHold still need their databases. **Experimental.**

- Structures only: proteins are folded with ESMFold first (≤ 400 aa; longer ones
  get no FoldSeek evidence).
- Hypothetical genes only, 4 at a time (`PHAGEFACTOR_WEB_WORKERS`); results cached
  by sequence in `<db_root>/web_cache` (shared, never pruned); failures are retried
  on the next run.
- Same scoring as local mode; hit domains from the taxid (taxdump).

It runs on the machine that launched the workflow (the login node on a cluster),
since compute nodes usually have no internet. On a workstation, prefer
`search_mode: local` with pdb100 + afdb-swissprot (~5.5 GB); `phagefactor check`
lists what is installed. `search_mode: none` curates on pharokka + phold only.

## Quick selector

| You have… | Input | Runner | Backend |
|---|---|---|---|
| Per-prophage genomes, cluster + DBs | genome | `--slurm` | `local` |
| Per-prophage genomes, workstation, small DBs | genome | local | `local` |
| Bulk CDS / pangenome proteins, cluster | protein / pangenome | `--slurm` | `local` |
| eggNOG-annotated source genome (citro) | protein (pre-CDS sub-mode) | either | `local` |
| No FoldSeek databases (Pharokka/PHold DBs only) | any | either | `webapi` |
