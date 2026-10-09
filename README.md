# phageFACTor

**Structure-aware functional annotation of prophage (phage) genomes.**
*FACT = Functional Annotation and Curation Tool.*

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

Validated on *Campylobacter*, *Citrobacter* and OMM12-consortium prophages.

---

## Where to run what

| | how | databases |
|---|---|---|
| **Cluster** (production) | `bash submit_all.sh` | afdb50 + pdb100 + swissprot |
| **Local** (testing, reviewers) | `bash run_local.sh` | pdb100 + swissprot (~5.5 GB) |

afdb50 needs ~60 GB and a large-memory node, so full-depth annotation is a
cluster job. `SEARCH_MODE=webapi` searches without local databases (it folds
queries with ESMFold first); the public ESMFold service is intermittent, so
treat it as experimental.

- **Genome mode**: one nucleotide FASTA per prophage → Pharokka → Phold → FoldSeek.
- **Protein mode**: bulk protein FASTA (all treated as hypothetical) → Phold proteins
  → FoldSeek. Supports a *pre-CDS-computed* sub-mode that rebuilds a GenBank from
  source-genome coordinates (for eggNOG-annotated inputs).
- **Pangenome sub-mode**: a bulk **nucleotide** CDS set / clustered pan-proteome
  (e.g. a Roary representative FASTA). Drop it in `input/fasta/`; it is translated
  (clean ids + id map), run as protein mode, and a **homology-reference set**
  (`pangenome_reference.{csv,xlsx,faa}`) is built for transferring annotations to a
  larger protein collection by homology (DIAMOND/MMseqs2). See [`docs/modes.md`](docs/modes.md).
- **Search mode** is chosen at step 02: local `foldseek` against installed DBs, or the
  public FoldSeek web server (ProstT5 runs server-side, no DB download).

See [`docs/modes.md`](docs/modes.md) for the full matrix.

---

## Install

```bash
# 1. get the code
git clone https://github.com/TMs-code/phagefactor.git
cd phagefactor

# 2. environment (micromamba recommended; conda works for non-org users)
micromamba env create -f environment.yml      # creates 'phagefactor'
micromamba activate phagefactor

# 3. databases (one-time) — see docs/databases.md
#    python scripts/tools/check_databases.py   # what you have / what is missing

# 4. configure for your machine — see "Configure" below
```

The repo ships path-free (blank database fields, generic SLURM settings). You
adapt it with **environment variables** — no need to touch the committed
`config/config.yaml` / `config.sh`. The analysis scripts never need editing.

---

## Configure (adapt to your machine)

**Recommended: a personal env file.** Put your machine's values in one small
`env_<name>.sh` (gitignored — keep paths/email out of the repo) and `source` it
before running, leaving the shipped config untouched. Full details:
[`docs/configure.md`](docs/configure.md).

```bash
# env_mymachine.sh  —  source before submit_all.sh
export PHAGEFACTOR_INPUT=/data/jobs/myjob/input      # holds fasta/ + prophage_list.txt
export PHAGEFACTOR_RUN_BASE=/data/jobs/myjob         # outputs -> <base>/run, run_002, ...
export PHAGEFACTOR_DB_ROOT=/data/databases           # -> <root>/pharokka_db, <root>/foldseek_dbs
export FOLDSEEK_DB_ROOT="${PHAGEFACTOR_DB_ROOT}/foldseek_dbs"
export SLURM_EMAIL=you@example.org                   # blank = no notifications
export MAMBA_ENV=phagefactor                         # env holding pharokka+phold+foldseek
export SCRATCH_TMPDIR=/local/scratch/tmp             # fast local disk for temp files
```
```bash
source env_mymachine.sh && bash submit_all.sh
```

Every setting resolves as **env var > `config` value > default**, so you can also
edit `config/config.yaml`/`config.sh` directly if you prefer.

| What | env var | `config` fallback | Default |
|---|---|---|---|
| Input folder | `PHAGEFACTOR_INPUT` | `paths.input_dir` | `<repo>/input` |
| Run-dir parent | `PHAGEFACTOR_RUN_BASE` | `paths.run_base` | `<repo>/runs` |
| FoldSeek DB root | `FOLDSEEK_DB_ROOT` | `databases.foldseek_db_root` | — |
| Pharokka / FoldSeek DB parent | `PHAGEFACTOR_DB_ROOT` | `config.sh DB_ROOT` | `<repo>/databases` |
| Email | `SLURM_EMAIL` | `slurm.email` | blank |
| Conda/mamba env | `MAMBA_ENV` | `config.sh` | `phagefactor` |
| Scratch dir | `SCRATCH_TMPDIR` | `config.sh` | `/tmp` |

Outputs go to a **versioned** run dir under `run_base` — `run/`, then `run_002/`,
… — a new submission **never overwrites** a previous run. To reuse one instead:
`RESUME=1 bash submit_all.sh` (most recent run) or `PHAGEFACTOR_RUN_DIR=/path` (a
specific run); add `SKIP_PHAROKKA=1 SKIP_PHOLD=1` to skip already-done steps.

**Cluster notes:** partition/QOS default to `common`/`fast`; phold is CPU by default (GPU needs a *verified* CUDA env or it
silently falls back); the local-FoldSeek step (02) loads the afdb50 index into RAM
and requests **256 GB** — on smaller nodes use WebAPI mode or the low-RAM option in
[`docs/databases.md`](docs/databases.md).

---

## Run on a cluster (SLURM) — the normal way

`submit_all.sh` (repo root) is the single driver. It auto-detects genome vs protein mode
from the input folder, submits the SLURM steps in `steps/`, and chains them with
dependencies:

```bash
# point at the code once (stays pristine); data + outputs live elsewhere:
cd /path/to/phagefactor
export PHAGEFACTOR_INPUT=/path/to/job/input        # holds fasta/ + prophage_list.txt
export PHAGEFACTOR_RUN_BASE=/path/to/job            # runs land in <base>/run, run_002, ...
bash submit_all.sh                                  # pharokka → phold → FoldSeek → compare → curate → output
```

`submit_all.sh` prints the resolved run dir and the optional follow-up to add phynteny +
C1/Cro synteny (already wired with `--chdir` so its logs join the run):

```bash
sbatch --chdir=<run_dir> --dependency=afterok:<jobid> steps/05_phynteny.sh
# review the few flagged genes:  <run_dir>/04_output/curation/review_suggested.csv
```

`SEARCH_MODE`: `local` (default) · `webapi` (no local DBs, experimental) ·
`none` (pharokka + phold only). `python scripts/tools/check_databases.py` shows
what is installed.

---

## Run on one machine (no SLURM) — testing

`run_local.sh` runs the same steps sequentially. Two *Campylobacter* prophages ship
in `example_data/genome_2campy/`:

```bash
bash steps/00b_setup_databases.sh --pharokka --pdb100 --swissprot   # ~7 GB, one-time
cp -r example_data/genome_2campy/input/* input/
FOLDSEEK_DBS="afdb-swissprot pdb100" bash run_local.sh
```

Without afdb50 the annotation rate is lower than a cluster run — this is for
checking the pipeline works, not for production. `bash example_data/run_smoke_test.sh`
does the same and verifies the output schema.

See [`example_data/README.md`](example_data/README.md) for all input modes and the
column reference.

---

## Outputs

`<run_dir>/04_output/final_annotations_table.{csv,xlsx}` — one row per CDS, with the
curated `final_product`/`short_name`, the **evidence category**, and full traceability
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
`05_phynteny/final_with_synteny.gb`. In **pangenome sub-mode** there is no GenBank;
instead `04_output/pangenome_reference.{csv,xlsx,faa}` is produced (one row per
group + a DIAMOND/MMseqs-ready FASTA).
Full column reference: [`docs/outputs.md`](docs/outputs.md).

---

## How annotations are decided

The curation engine is deliberately conservative: ambiguous genes are routed to
review rather than auto-annotated. Core principles:

- **A label is not a function.** A FoldSeek-only name that is only a domain/family
  label ("SHOCT domain-containing protein") is not used as the final annotation;
  it stays in `fs_description` / `fs_top3`.
- **Reference sub-databases are weighed by identity.** A Phold hit to NetFlaX, VFDB
  or CARD below 40% identity yields to a CONFIDENT FoldSeek call that names a function.
- **Structure titles are not names.** PDB titles are reduced to the protein name
  they contain, or replaced by a top-3 hit that names a function.
- **Stronger evidence wins.** A name repeated in the FoldSeek top-3 replaces the
  best hit only if its alignment score is at least as high.

All vocabulary lives in `scripts/lib/lexicon.py` (word-lists, order irrelevant) and
`scripts/lib/rules.py` (ordered rules: name upgrades, PHROG category, `short_name`).
Thresholds are in `config/config.yaml` ([`docs/configure.md`](docs/configure.md)).
Details: [`docs/curation_gates.md`](docs/curation_gates.md),
[`docs/code_lists.md`](docs/code_lists.md), diagrams in
[`docs/architecture.md`](docs/architecture.md).

---

## Citing

If you use phageFACTor, please cite it (see [`CITATION.cff`](CITATION.cff)) **and** its
dependencies: Pharokka, Phold, FoldSeek, Phynteny, and (for pre-computed inputs)
eggNOG-mapper.

## License

MIT — see [`LICENSE`](LICENSE).
