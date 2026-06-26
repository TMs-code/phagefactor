# Modes: genome vs protein × local vs WebAPI

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
  window coordinates (+ a rich gene-metadata CSV with absolute coords), step 05
  rebuilds a GenBank so genome-style outputs work even though there was no
  Pharokka run (used for the citro phiNP/phiSM eggNOG inputs).

## Search mode (chosen at step 02)

### Local (`02d_foldseek_3di.py`) — default, cluster-scale
Searches the Phold 3Di tokens against installed DBs (PDB100, AFDB-SwissProt,
AFDB50/BakTFold). Builds the query DB via the official ProstT5 `tsv2db` recipe
(no `_ca` file). Per-DB taxonomy + taxon filtering, per-residue masking. Needs the
DBs on disk — see [`databases.md`](databases.md).

### WebAPI (`02w_foldseek_webapi.py`) — laptop-friendly, no local DB
Submits the AA FASTA to `search.foldseek.com`; ProstT5 + search run server-side.
Ticket/poll with resubmit-on-timeout. Output: `webapi_hits.m8`.

> **Known gap (v1 to close):** the WebAPI m8 is not yet auto-fed through
> `foldseek_scoring` to produce `best_hit.csv`/`top3.csv`. Until a `--to-best-hit`
> bridge is added, the WebAPI path needs a manual parse/score step before `03`.

## Quick selector

| You have… | Input mode | Search mode |
|---|---|---|
| Per-prophage genomes, cluster + DBs | genome | local |
| Per-prophage genomes, laptop, no DBs | genome | WebAPI |
| Bulk CDS / pangenome proteins, cluster | protein | local |
| eggNOG-annotated source genome (citro) | protein (pre-CDS sub-mode) | local or WebAPI |
