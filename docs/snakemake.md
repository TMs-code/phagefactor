# Running phageFACTor

phageFACTor is a Snakemake workflow (`workflow/Snakefile`) launched with
`phagefactor run`. Snakemake works out from the files themselves what still has
to run:

- **New analysis:** everything runs.
- **Edited curation rule, lexicon entry or threshold:** only the steps that use it re-run.
- **Crashed or cancelled run:** relaunch the same command and it resumes where it stopped.
- **Same machine and cluster:** one command for both (`--slurm`).

---

## 1. Install (once)

```bash
git clone https://github.com/TMs-code/phagefactor.git
export PATH="$PWD/phagefactor/bin:$PATH"            # add to ~/.bashrc
micromamba create -y -f phagefactor/envs/pf_snakemake.yaml   # snakemake >= 8.20 + SLURM plugin
micromamba activate pf_snakemake
phagefactor install-envs                             # envs 'phagefactor' (tools) and 'phynteny'
phagefactor databases /path/to/db_dir --all          # see databases.md (afdb50: cluster)
```

Snakemake runs in its own small environment (the *driver*); each job runs in the
tool environment with `micromamba run -n <env>` (or `conda run`). With
`tools_env: auto` (default), jobs use the active environment when it contains
pharokka, phold and foldseek (a conda install of phageFACTor), and the env named
`phagefactor` otherwise; likewise `phynteny_env` and `phynteny`. If micromamba is
not found, set its path with `micromamba:` in `config/workflow.yaml`.

---

## 2. Run

```bash
phagefactor init /path/Jobs/MyProject       # once: config/ + input/
# input/fasta/ (+ input/prophage_list.txt in genome mode); edit config/*.yaml

tmux new -s pf                              # the driver must outlive your login session
micromamba activate pf_snakemake            # NOT the tool env: the jobs switch to it themselves
cd /path/Jobs/MyProject
phagefactor run --slurm -n                  # dry run: list what WOULD run, and why
phagefactor run --slurm                     # submit
# detach: ctrl-b d      come back: tmux attach -t pf

phagefactor run                             # or: everything on this machine
```

`phagefactor run` is `snakemake --snakefile <code>/workflow/Snakefile --profile
<code>/profiles/{local,slurm}`; every other argument goes to Snakemake. With
Snakemake 8, add `-r` to the dry run to see the reasons.

On SLURM, keep Snakemake running on the login node in `tmux` or `screen`. It
only submits jobs and waits for them; it uses no CPU. If the session dies, jobs
already submitted keep running but nothing new is submitted; relaunch the same
command to continue (jobs that were still running are redone).

### Safeguards

- **Launch check** (also with `-n`): the tools and databases of the steps still to
  run must be found, else the run stops with what is missing. Skip:
  `--config preflight=false`.
- **Early abort:** after 5 failed Pharokka / Phold jobs and no success, the rest
  of that step exits at once.
- **Mail** (`slurm_email`): on failure of single steps, not of each Pharokka /
  Phold job.

### Settings

`config/workflow.yaml`: what to run and where it goes. `config/config.yaml`: the
science (thresholds, databases, curation). Machine settings (env names, SLURM
account, node features) go in `~/.config/phagefactor/settings.yaml`. For each key,
the value is taken from `--config key=value`, else the project file, else your
machine settings, else the environment variable, else the default. Reference:
[`configure.md`](configure.md).

Each step's `runtime` must fit within the MaxWall of its QOS. Check with
`sacctmgr show qos format=Name,MaxWall`. `phold` / `phold_proteins` ask for
120 min under `fast`, exactly the 2 h ceiling.

Examples:

```bash
phagefactor run --slurm --config run_dir=runs/webapi search_mode=webapi
phagefactor run --slurm --config phynteny=false     # skip phynteny
phagefactor run --slurm --config 'resources={"foldseek_search":{"runtime":720}}'
```

### Interactive report

Every run ends with `<run_dir>/04_output/report.html`: one self-contained file (no
internet needed) with a genome map per prophage coloured by function, a filterable
table and the evidence behind each call. It is built on the machine running
Snakemake (a few seconds, no SLURM job), after phynteny when phynteny runs.

### Phynteny runs by default

Phynteny + synteny integration runs after curation
unless `phynteny: false`. It is flagged so you always know:

- at the start: `[phageFACTor] PHYNTENY: ON (default) -> <run>/05_phynteny/final_annotations_integrated.csv`
- at the end: `PHYNTENY RAN: ...`, or `PHYNTENY WAS SKIPPED: <reason>` when no
  GenBank could be built (e.g. plain protein mode).

`04_output/final_annotations_table.csv` is unchanged; the phynteny version is the
separate `05_phynteny/final_annotations_integrated.csv`. Phynteny's ESM2 model
(~2.5 GB) is downloaded once into `phynteny_cache` (default
`~/.cache/phagefactor/phynteny`) and shared by all runs.

### One run dir per analysis

Relaunching on a run dir redoes only what changed. For a separate analysis, give
a new `run_dir`; two run dirs can run at the same time from the same project.

---

## 3. What re-runs when

Snakemake compares what each output was made from with what it would be made
from now:

| You change | Re-runs |
|---|---|
| curation code (`04_curate_annotations.py`, `05_build_output.py`) or `curate:` thresholds | curate (+ build output) |
| `lib/rules.py` (naming rules, also used by the FoldSeek vote) | FoldSeek post-processing → compare → curate |
| `lib/lexicon.py`, `lib/utils.py`, `lib/taxonomy.py`, installing the NCBI taxdump | FoldSeek post-processing → compare → curate |
| `foldseek:` scoring settings (tiers, `low_info_malus`, `host_boost`…), `lib/foldseek_scoring.py`, `02d_foldseek_3di.py` | FoldSeek post-processing → compare → curate |
| `compare:` thresholds, `03_compare_annotations.py` | compare → curate |
| `evalue_max`, the DB list, `prostt5_mask_threshold` | **FoldSeek search** → onward |
| one genome FASTA | that prophage's pharokka + phold → merge → FoldSeek → onward |
| one protein in a protein FASTA | the one phold batch that contains it → merge → onward |
| `input/overrides.tsv`, or the precomputed sub-mode inputs (`genome/source_genome.fasta`, `fasta/prophage_windows.csv`, `gene_metadata_rich.csv`) | curate (+ build output) |
| `foldseek_taxon_filter` (`databases:` in `config.yaml`) | FoldSeek post-processing → compare → curate |
| `scripts/lib/config.py` (shared paths + defaults) | FoldSeek post-processing → compare → curate → phynteny |
| `06_phynteny.py`, `07_integrate.py`, `lib/synteny_hint.py` | phynteny only |
| `09_report.py`, `lib/report_template.html` | the report only |

Phynteny reads the final table, so every row that re-runs curate also re-runs
phynteny, and the report is rebuilt after either. Editing only a step's shell
command in the Snakefile does not re-run finished work (profiles:
`rerun-triggers: [mtime, params, input, software-env]`).

The FoldSeek search (hours, 256 GB) re-runs only when its own inputs or search
settings change, never because the scoring or curation code changed. To rebuild
a deleted intermediate file, ask for its step explicitly (`-R compare`).

### Useful commands

```bash
phagefactor run --slurm -n               # dry run with the reason for each job (8.x: add -r)
phagefactor run --slurm -R curate        # force a step (and everything after it)
phagefactor run --slurm --until compare  # stop after a step
phagefactor run --slurm --unlock         # after a killed run: "Directory cannot be locked"
phagefactor run --dag | dot -Tsvg > dag.svg   # picture of the jobs
```

- **Step logs:** `<run_dir>/logs/smk/<step>.log`.
- **SLURM job logs:** `<run_dir>/.snakemake/slurm_logs/`.

### Plain Snakemake, as a module (since 1.2)

The repo follows the
[workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog) layout,
so it also runs without the `phagefactor` command, and other workflows can include it:

```bash
snakedeploy deploy-workflow https://github.com/TMs-code/phagefactor myjob --tag v1.2.0
cd myjob && snakemake --cores all --sdm conda
```

`snakedeploy` writes a `workflow/Snakefile` that loads this one from GitHub
(`module phagefactor: snakefile: github(...)`) and copies `config/`. The project
dir is the directory Snakemake runs in, as with `phagefactor run`. The scripts
are fetched into `.snakemake/phagefactor/code/`. A host workflow passes its
settings with `config:`: the `config.yaml` sections (`project`, `foldseek`,
`curate`, ...) are merged over `config/config.yaml`, and every other key is a
`workflow.yaml` setting. Details: [`config/README.md`](../config/README.md).

`--sdm conda` builds `envs/phagefactor.yaml` and `envs/phynteny.yaml` and runs
each job in its rule's env. Without it, jobs use `tools_env` / `phynteny_env` as
above. Keep one route per run dir: the profiles' `software-env` rerun trigger can
treat a switch as a software change and re-run finished steps.

---

## 4. Run dirs made by version 1.0

They are continued as they are (`--config run_dir=/path/to/run`): the dry run
lists FoldSeek post-processing → compare → curate → phynteny → report, and keeps
pharokka, phold and the search. If it also lists `pharokka`, `split_proteins` or
the search although the inputs did not change (input files copied again after the
run), mark those steps as done once, then dry-run again:

```bash
phagefactor run --slurm --config run_dir=/path/to/run --touch --until foldseek_search
```

---

## 5. Tests

```bash
micromamba activate pf_snakemake          # has snakemake + pandas/biopython for the mocks
bash tests/run_workflow_tests.sh           # ~2-3 min, needs no DB / GPU / SLURM
```

The tests run every mode with mock tools (`tests/mock_bin/`). They check:

- every input and search mode;
- the re-run rules in the table above;
- failure + resume, and stale batch handling;
- tool environment lookup;
- `phagefactor init` / `run` with a read-only code dir (nothing written to it);
- taking over a run dir without Snakemake records;
- use as a module loaded from a git source (snakedeploy layout), `--sdm conda`
  (dry run, when conda is on PATH) and the catalog's rule graph from `.test/`.

Run them after any change to the Snakefile or the scripts.
