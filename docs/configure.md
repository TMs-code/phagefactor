# Configuring phageFACTor

The code folder (a git clone, or a read-only conda install) holds no data and is
never written to. Everything about an analysis lives in a **project dir**:

```
myjob/
  config/config.yaml      paths, databases, thresholds (science)
  config/workflow.yaml    what to run and where (modes, envs, SLURM, resources)
  input/                  fasta/, prophage_list.txt, overrides.tsv
  runs/run/               outputs + logs (created by the run)
```

`phagefactor init myjob` creates it from the defaults shipped with the code; a git
clone of the repo is itself a project dir. Launch every run from the project dir
(`cd myjob && phagefactor run ...`). A file missing from `config/` falls back to the
shipped default.

Resolution order for every setting:
**`--config key=value` > project file > your machine settings > environment variable > default.**

## Your machine settings (`~/.config/phagefactor/settings.yaml`)

Settings that belong to a machine or a cluster account, not to a project or to the
code: the database folder, the conda env names, SLURM partition/QOS/account,
email, node features. They go in one personal file that every project reads, so
neither the code nor the projects carry them. It takes the same keys as
`config/workflow.yaml`, plus `db_root` (written by `phagefactor databases <dir>`):

```yaml
# ~/.config/phagefactor/settings.yaml
db_root: /path/to/phagefactor_db
tools_env: my_tools_env
slurm_partition: common
slurm_qos: fast
slurm_email: me@example.org
resources:                       # per step; only the keys given change
  foldseek_search: {constraint: avx512}
```

`constraint` asks SLURM for nodes with those features (`sbatch --constraint`).
FoldSeek's pre-filter gives slightly different hits on different CPU types, so
runs meant to be compared should search on the same node type.

## Databases (once per machine)

```bash
phagefactor databases /data/phagefactor_db --all     # or a subset of flags
```

records `db_root: /data/phagefactor_db` in `~/.config/phagefactor/settings.yaml`, so
every project finds `pharokka_db/`, `phold_db/` and `foldseek_dbs/` there. A project
can point elsewhere with `databases.db_root` (or `pharokka_db`, `phold_db`,
`foldseek_db_root`, `taxdump`) in its `config.yaml`. Without any of these, `<project>/databases`
is used. `phagefactor where` prints the resolved paths, `phagefactor check` what is
installed. When `phold_db/` does not exist, phold uses its own default location.

## Run settings (`config/workflow.yaml`)

| Key | Meaning | Default |
|---|---|---|
| `input_dir` / `run_dir` | input folder / run dir (relative to the project dir) | `paths.*` of config.yaml: `input`, `runs/run` |
| `input_mode` | `auto`, `genome`, `protein`, `pangenome` | `auto` |
| `search_mode` | `local`, `webapi`, `none` | `local` |
| `foldseek_dbs` | subset of the local DBs, e.g. `"afdb-swissprot pdb100"` | all |
| `phynteny` | phynteny + synteny integration | `true` |
| `tools_env` / `phynteny_env` | conda env of the tools; `auto` = the active env if it has them, else `phagefactor` / `phynteny` | `auto` |
| `slurm_partition`, `slurm_qos`, `slurm_account`, `slurm_email`, ... | SLURM submission | `common`, `fast`, none, none |
| `resources:` | CPUs / memory / time per step | see the file |
| `preflight` | check tools and databases of the pending steps before submitting | `true` |

Every key also has an environment variable (listed in the file), so an existing
`env_<project>.sh` keeps working.

## Annotation settings (`config/config.yaml`)

These change annotations: change them deliberately. Editing one re-runs only the
steps that read it.

| Key | Default | Meaning |
|---|---|---|
| `project.host_genus` | blank | host genus for the same-host flag and the host-derived-moron gate (env `PHAGEFACTOR_HOST_GENUS`) |
| `foldseek.confident_evalue` / `confident_score` / `confident_top3_score` | 1e-3 / 200 / 300 | CONFIDENT tier |
| `foldseek.good_evalue` / `good_score` / `good_top3_evalue` | 0.01 / 200 / 0.05 | GOOD tier |
| `foldseek.borderline_evalue` / `borderline_score` / `borderline_scoreonly_evalue` | 0.1 / 90 / 0.5 | BORDERLINE tier (below = WEAK) |
| `foldseek.low_info_malus` | 0.60 | ranking malus for domain/family/repeat labels and gene labels |
| `foldseek.host_boost` | 1.20 | same-host ranking boost (env `PHAGEFACTOR_HOST_BOOST`) |
| `curate.subdb_close_min_fident` | 0.4 | Phold identity at which a NetFlaX/VFDB/CARD hit is "close" (see `curation_gates.md`) |

## The run dir

One run dir per analysis. Relaunching on the same run dir redoes only what changed
(an input FASTA, a rule, a threshold) or what failed; nothing is overwritten
otherwise. For a second analysis of the same project use another run dir:
`phagefactor run --config run_dir=runs/strict_evalue`. Each run dir contains
`00c_pharokka/ ... 05_phynteny/` and `logs/`; the Snakemake job logs are in
`logs/smk/`. Re-run rules: [`snakemake.md`](snakemake.md).

## Environment variables

| Variable | Sets |
|---|---|
| `PHAGEFACTOR_INPUT` / `PHAGEFACTOR_RUN_DIR` | input folder / run dir |
| `PHAGEFACTOR_DB_ROOT` | parent of `pharokka_db`, `phold_db`, `foldseek_dbs` |
| `PHAROKKA_DB` / `PHOLD_DB` / `FOLDSEEK_DB_ROOT` / `PHAGEFACTOR_TAXDUMP` | one database location |
| `PHAGEFACTOR_CONFIG` | a settings file other than `<project>/config/config.yaml` |
| `SCRATCH_TMPDIR` | fast local disk for temporary files |
| `MAMBA_ENV` / `PHYNTENY_ENV` | `tools_env` / `phynteny_env` |
| `PHAGEFACTOR_WEB_CACHE` | web-API cache folder (default `<db_root>/web_cache`) |
| `PHAGEFACTOR_WEB_WORKERS` | proteins searched at once in web-API mode (default 4) |
