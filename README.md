# phageFACTor

**Structure-aware functional annotation of prophage (phage) genomes.**
*FACT = Filter / Assign / Curate Tool.*

phageFACTor rescues the "hypothetical protein" fraction of prophage annotations by
combining three evidence layers and reconciling them with an auditable curation
engine:

1. **Pharokka** — gene calling (Phanotate) + PHROG / specialised-DB annotation.
2. **Phold** — ProstT5 3Di structural tokens + structural search vs phage DBs.
3. **Custom FoldSeek** — the 3Di tokens searched against general structural DBs
   (PDB100, AFDB-SwissProt, AFDB50) for hits Phold misses.

A comparison + curation step (**the FAC core**) reconciles Phold vs custom-FoldSeek
per gene, assigns an evidence category, and emits a curated table + updated GenBank.
Optional: **Phynteny** synteny categories and a positional C1/Cro switch heuristic.

> Status: research code, first repo-grade draft. The scientific logic is mature and
> validated on *Campylobacter*, *Citrobacter* and *OMM12 consortia* prophages.

---

## The 2 × 2 it supports

|                 | **Local FoldSeek** (DBs on disk) | **WebAPI FoldSeek** (no local DB) |
|-----------------|----------------------------------|-----------------------------------|
| **Genome mode** (nucleotide FASTAs) | full pipeline, cluster-scale | laptop-friendly, no big DB |
| **Protein mode** (bulk CDS / pre-defined GBK) | full pipeline | laptop-friendly |

- **Genome mode**: one nucleotide FASTA per prophage → Pharokka → Phold → FoldSeek.
- **Protein mode**: bulk protein FASTA (all treated as hypothetical) → Phold proteins
  → FoldSeek. Supports a *pre-CDS-computed* sub-mode that rebuilds a GenBank from
  source-genome coordinates (for eggNOG-annotated inputs such as citro).
- **Search mode** is chosen at step 02: local `foldseek` against installed DBs, or the
  public FoldSeek web server (ProstT5 runs server-side, no DB download).

See [`docs/modes.md`](docs/modes.md) for the full matrix.

---

## Install

```bash
# 1. environment (micromamba recommended; conda works for non-org users)
micromamba env create -f environment.yml      # creates 'phagefactor'
micromamba activate phagefactor

# 2. databases (one-time) — see docs/databases.md
#    pharokka_db, and FoldSeek DBs (skip if using WebAPI search mode)

# 3. point config at your machine
#    edit config/config.yaml -> databases:  (or export FOLDSEEK_DB_ROOT)
#    everything else resolves relative to the repo, or to $PHAGEFACTOR_ROOT
```

No personal paths: paths resolve from `$PHAGEFACTOR_ROOT` (or the repo
root) and the single `config/config.yaml`. See [`docs/databases.md`](docs/databases.md)
for the FoldSeek DB layout and the **no-index / low-RAM** option.

---

## Quick start (shipped example, WebAPI search — no local DB needed)

Two *Campylobacter* prophages ship in `example_data/genome_2campy/`. After installing
the env (above), the smoke test stages them and runs the pipeline end-to-end against the
public FoldSeek web server (no database download), then checks the output:

```bash
cd phagefactor
bash example_data/run_smoke_test.sh
# -> 04_output/final_annotations_table.csv  (+ 04_output/curation/review_suggested.csv)
```

Or run it yourself in WebAPI mode:

```bash
cp -r example_data/genome_2campy/input/* input/
SEARCH_MODE=webapi bash submit_all.sh
```

See [`example_data/README.md`](example_data/README.md) for all input modes and the
column reference.

---

## Run on a cluster (SLURM, local FoldSeek)

`submit_all.sh` (repo root) is the single driver. It auto-detects genome vs protein mode
from `input/`, submits the SLURM steps in `steps/`, and chains them with dependencies:

```bash
# genome mode: per-prophage FASTAs in input/fasta/, names in input/prophage_list.txt
bash submit_all.sh                       # pharokka → phold → FoldSeek → compare → curate → output
# the run prints the optional follow-up to add phynteny + C1/Cro synteny:
sbatch steps/05_phynteny.sh
# review the few flagged genes:
#   04_output/curation/review_suggested.csv
```

Set `SEARCH_MODE=webapi` to skip local FoldSeek DBs. Snakemake is a next planned step.

---

## Outputs

`04_output/final_annotations_table.{csv,xlsx}` — one row per CDS, with the curated
`final_product`/`short_name`, the **evidence category**, and full traceability
(Pharokka / Phold / FoldSeek columns side by side). Evidence categories:

| Category | Meaning |
|---|---|
| **both agree** | Phold and custom-FoldSeek concur |
| **merged** | combined into one more-specific call |
| **structural-only** | custom-FoldSeek structural homology, no Phold corroboration |
| **phold-only** | Phold call, no FoldSeek support |
| **still-hypothetical** | neither found informative evidence |

Plus `updated_prophages.gb` (GenBank with `/product`, `/function`, `/note` provenance),
and after the phynteny step, `05_phynteny/final_annotations_integrated.csv` +
`05_phynteny/final_with_synteny.gb`.
Full column reference: [`docs/outputs.md`](docs/outputs.md).

---

## How annotations are decided

The curation engine (Rules A–E, re-routing gates, promiscuous-fold and eukaryotic
filters, phage-boost factor) is documented in
[`docs/curation_gates.md`](docs/curation_gates.md). It is deliberately conservative:
ambiguous genes are routed to `needs_review` rather than auto-annotated.

---

## Citing

If you use phageFACTor, please cite it (see [`CITATION.cff`](CITATION.cff)) **and** its
dependencies: Pharokka, Phold, FoldSeek, Phynteny, and (for pre-computed inputs)
eggNOG-mapper.

## License

MIT — see [`LICENSE`](LICENSE).
