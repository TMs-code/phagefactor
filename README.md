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
per gene, assigns an evidence category, and emits a curated table, an updated GenBank
and an **interactive HTML report** (genome maps coloured by function, and for every
gene the evidence behind its call). Optional: **Phynteny** synteny categories and a
positional C1/Cro switch heuristic.

Developed on the prophages of the OMM12 bacterial consortium and the Philympics 2021
prophage benchmark set ([Roach et al. 2021](https://f1000research.com/articles/10-758/v2)),
with manual curation of the resulting annotations.

---

## How it runs

phageFACTor is a Snakemake workflow driven by one command, `phagefactor`. It
re-runs only what changed: edit a curation rule and only curation re-runs; a
crashed run resumes where it stopped.

| | command | databases |
|---|---|---|
| **Cluster** (production) | `phagefactor run --slurm` | afdb50 + pdb100 + swissprot |
| **One machine** (testing, small datasets) | `phagefactor run` | pdb100 + swissprot (~5.5 GB) |
| **No FoldSeek databases** | `--config search_mode=webapi` | none (public FoldSeek server) |

afdb50 (~60 GB, plus a ~350 GB index loaded into RAM on a 256 GB node) makes
full-depth annotation a cluster job; its `_ca` files can be deleted.
`search_mode=webapi` searches the same three databases on the public FoldSeek
server, after folding with ESMFold (proteins ≤ 400 aa).

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
  public FoldSeek web server (no FoldSeek DB download).

See [`docs/modes.md`](docs/modes.md) for the full matrix.

---

## No cluster or databases? Request a run

I can run phageFACTor for you on our HPC, with the full databases and the afdb50
index, at no cost. Open an
[annotation request](https://github.com/TMs-code/phagefactor/issues/new?template=annotation_request.yml)
with the number of genomes (or proteins) and the input mode. **Do not attach
sequences to the issue**: I will reply with how to send them, plus the host genome if
you have one. Best effort, as cluster time allows.

---

## Install

Three ways in. All of them use the same workflow and conda-forge + bioconda
packages only (no Anaconda `defaults` channel). conda, mamba and micromamba all work.

**Bioconda** (the `phagefactor` command with the tools in the same env):

```bash
conda create -n phagefactor -c conda-forge -c bioconda phagefactor
conda activate phagefactor
phagefactor install-envs --phynteny   # phynteny has its own env
phagefactor databases /data/phagefactor_db --pharokka --phold --pdb100 --swissprot --taxdump
```

**Snakemake workflow** ([workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog)
layout, to run with plain Snakemake or combine with other workflows):

```bash
snakedeploy deploy-workflow https://github.com/TMs-code/phagefactor myjob --tag v1.2.0
cd myjob     # edit config/config.yaml (databases.db_root), put the FASTAs in input/fasta/
snakemake --cores all --sdm conda
```

Snakemake builds the tool envs itself. See [`config/README.md`](config/README.md).

**Git clone** (development):

```bash
git clone https://github.com/TMs-code/phagefactor.git
export PATH="$PWD/phagefactor/bin:$PATH"           # the `phagefactor` command

micromamba env create -f phagefactor/envs/pf_snakemake.yaml   # snakemake (driver)
micromamba activate pf_snakemake
phagefactor install-envs        # tool env (pharokka, phold, foldseek) + phynteny env
phagefactor databases /data/phagefactor_db --pharokka --phold --pdb100 --swissprot --taxdump
phagefactor check               # what is installed / missing
```

`phagefactor databases <dir>` records `<dir>` in `~/.config/phagefactor/settings.yaml`,
so every project uses the same databases. afdb50 and its index are cluster
downloads: [`docs/databases.md`](docs/databases.md).

---

## Run

A **project dir** holds the settings, the input and the runs.

```bash
phagefactor init myjob           # myjob/config/{config,workflow}.yaml + myjob/input/
cp my_prophages/*.fasta myjob/input/fasta/
ls myjob/input/fasta | sed 's/\.fasta$//' > myjob/input/prophage_list.txt   # genome mode
cd myjob
phagefactor run -n               # dry run: what would run
phagefactor run --slurm          # SLURM, from the login node (inside tmux/screen)
phagefactor run                  # or on this machine
```

Results go to `myjob/runs/run/`. Relaunching the same command after an edit or
a crash redoes only what is needed; a new analysis is a new project dir (or
`--config run_dir=runs/<name>`). Set the host genus, partition, QOS and email in
`config/config.yaml` and `config/workflow.yaml`; any key can also be given on the
command line (`--config search_mode=none`). Details:
[`docs/configure.md`](docs/configure.md), [`docs/snakemake.md`](docs/snakemake.md).

phold runs on CPU by default (GPU needs a *verified* CUDA env: `phold_gpu: true`).
Without a 256 GB node, search pdb100 + swissprot only, or use `search_mode=webapi`
([`docs/databases.md`](docs/databases.md)).

**Try it:** `bash example_data/run_smoke_test.sh` runs the two *Campylobacter*
prophages of `example_data/genome_2campy/` with the small databases and checks the
output. Without afdb50 the annotation rate is lower than a cluster run: this checks
that the pipeline works, it is not a production run. All input modes and the column
reference: [`example_data/README.md`](example_data/README.md).

---

## Outputs

`<run_dir>/04_output/final_annotations_table.{csv,xlsx}` — one row per CDS, with the
curated `final_product`/`short_name`, the **evidence category**, and full traceability
(Pharokka / Phold / FoldSeek columns side by side). Evidence categories:

| `annotation_source` | Meaning |
|---|---|
| `both agree` | Phold and FoldSeek name the same function |
| `merged` | the two calls combined into one more specific name |
| `foldseek` | FoldSeek-only: named from the FoldSeek hits, Phold had no informative call |
| `phold` | Phold-only: named from Phold, no FoldSeek support |
| `no_informative_hit` / `no_hit` | still hypothetical: hits without a function, or none |

`<run_dir>/04_output/report.html` — the interactive report, one self-contained file that
opens in any browser offline: a genome map per prophage coloured by function, a
searchable/filterable table, and per gene the evidence behind the call (Phold, FoldSeek
top-3, curation reason, phynteny). For a run made before version 1.1:
`python scripts/09_report.py --run_dir <run_dir>`.

Plus `updated_prophages.gb` (GenBank with `/product`, `/function`, `/note` provenance),
and from the phynteny step, `05_phynteny/final_annotations_integrated.csv` +
`05_phynteny/final_with_synteny.gb`. In **pangenome sub-mode** there is no GenBank;
instead `04_output/pangenome_reference.{csv,xlsx,faa}` is produced (one row per
group + a DIAMOND/MMseqs-ready FASTA).
Full column reference: [`docs/outputs.md`](docs/outputs.md).

---

## How annotations are decided

Every call is made by explicit rules and is traceable: the final table records the
rule that decided each gene (`curation_action`) and why (`curation_explanation`).
Ambiguous genes go to review instead of being guessed.

- **Evidence is weighed, not just collected.** Each FoldSeek hit is rated
  CONFIDENT / GOOD / BORDERLINE / WEAK from its e-value, score and agreement among
  the top-3 hits. A BORDERLINE hit never names a gene on its own.
- **Phage-aware ranking.** Hits naming a phage function are ranked up (×2 / ×1.5),
  uninformative labels down; the boost never lifts a weak hit above a confident one.
- **Taxonomy-aware.** Each hit's domain comes from the FoldSeek lineage or NCBI
  taxid. Eukaryotic hits, mostly look-alike folds, get no phage boost, never outrank
  a confident bacterial hit, and are demoted unless a bacterial hit names the same
  function.
- **Specific beats generic.** A generic Phold call
  ("transcriptional regulator") gives way to a confident specific FoldSeek name;
  complementary calls are merged ("tail protein (Late control protein)").
- **Conflicts are flagged, not resolved by guesswork.** Confident calls that
  disagree, or that place a protein on the head versus the tail, go to review.
- **A label is not a function.** Domain/family labels and structure titles are not
  used as names, except characterized domains (`lexicon.RELEVANT_DOMAINS`, sent to
  review) and labels renamed by a naming rule (`rules.UPGRADE_RULES`).
- **Reference sub-databases are weighed by identity.** A Phold hit to NetFlaX, VFDB
  or CARD below 40% identity yields to a confident FoldSeek call that names a function.

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
