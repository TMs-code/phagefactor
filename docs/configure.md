# Configuring phageFACTor for your machine

phageFACTor ships **path-free**: `config/config.yaml` and `config.sh` have blank
databases and generic SLURM defaults. You adapt it entirely with **environment
variables** — the committed config files never need editing.

Resolution order for every setting: **environment variable > value in `config` > built-in default.**

## The env file (recommended)

Create one small file per machine/project, keep it out of git (`env_*.sh` is
gitignored), and `source` it before submitting:

```bash
# env_mymachine.sh
export PHAGEFACTOR_INPUT=/data/jobs/myjob/input   # folder with fasta/ + prophage_list.txt
export PHAGEFACTOR_RUN_BASE=/data/jobs/myjob      # outputs -> <base>/run, run_002, ...
export PHAGEFACTOR_DB_ROOT=/data/databases        # -> <root>/pharokka_db and <root>/foldseek_dbs
export FOLDSEEK_DB_ROOT="${PHAGEFACTOR_DB_ROOT}/foldseek_dbs"
export SLURM_EMAIL=you@example.org                # blank = no notifications
export MAMBA_ENV=phagefactor                      # conda/mamba env with pharokka+phold+foldseek
export SCRATCH_TMPDIR=/local/scratch/tmp          # fast local disk for temp files
```

```bash
source env_mymachine.sh
bash submit_all.sh
```

Why an env file instead of editing `config.yaml`/`config.sh`: it keeps your
personal paths/email out of the repo, and the shipped config files stay identical
to what everyone else runs.

## Full variable reference

| Variable | Sets | Config fallback | Default |
|---|---|---|---|
| `PHAGEFACTOR_INPUT` | input folder (`fasta/` + `prophage_list.txt`) | `paths.input_dir` | `<repo>/input` |
| `PHAGEFACTOR_RUN_BASE` | parent of the versioned run dirs | `paths.run_base` | `<repo>/runs` |
| `PHAGEFACTOR_RUN_DIR` | exact run dir (pins/reuses one) | — | auto-versioned |
| `PHAGEFACTOR_DB_ROOT` | parent of `pharokka_db` + `foldseek_dbs` (bash) | `config.sh DB_ROOT` | `<repo>/databases` |
| `FOLDSEEK_DB_ROOT` | FoldSeek DB root (Python step 02) | `databases.foldseek_db_root` | — |
| `SLURM_EMAIL` | notification email | `slurm.email` | blank |
| `MAMBA_ENV` | conda/mamba env name | `config.sh MAMBA_ENV` | `phagefactor` |
| `SCRATCH_TMPDIR` | scratch dir for temp files | `config.sh` | `/tmp` |
| `PARTITION` / `QOS` | SLURM partition / QOS | `config.sh` | `common` / `fast` |

WebAPI search mode (`SEARCH_MODE=webapi`) needs no FoldSeek DB at all — skip the
DB variables.

## Annotation settings (scientific, in `config/config.yaml`)

These are not machine-specific, so they live in `config.yaml` rather than in an env
file. Change them only deliberately: they change annotations.

| Key | Default | Meaning |
|---|---|---|
| `foldseek.confident_evalue` / `confident_score` / `confident_top3_score` | 1e-3 / 200 / 300 | CONFIDENT tier |
| `foldseek.good_evalue` / `good_score` / `good_top3_evalue` | 0.01 / 200 / 0.05 | GOOD tier |
| `foldseek.borderline_evalue` / `borderline_score` / `borderline_scoreonly_evalue` | 0.1 / 90 / 0.5 | BORDERLINE tier (below = WEAK) |
| `foldseek.low_info_malus` | 0.60 | ranking malus for domain/family/repeat labels and gene labels |
| `foldseek.host_boost` | 1.20 | same-host ranking boost (env `PHAGEFACTOR_HOST_BOOST`) |
| `curate.subdb_close_min_fident` | 0.4 | Phold identity at which a NetFlaX/VFDB/CARD hit is "close" (see `docs/curation_gates.md`) |

## The run dir (outputs + logs)

Every submission writes to a **versioned** dir under `PHAGEFACTOR_RUN_BASE`:
`run/`, then `run_002/`, `run_003/`, … so a run **never overwrites** a previous
one. `submit_all.sh` picks the next free dir, creates it, exports
`PHAGEFACTOR_RUN_DIR`, and writes `<run_base>/.current_run` (a pointer that later
manual steps follow). Each run dir contains `00c_pharokka/ … 05_phynteny/` and
`logs/`.

## Resuming / starting from a step

| Goal | How |
|---|---|
| New run (default) | `bash submit_all.sh` → fresh `run_00N/` |
| Reuse the most recent run | `RESUME=1 bash submit_all.sh` |
| Reuse a specific run | `PHAGEFACTOR_RUN_DIR=/…/run_002 bash submit_all.sh` |
| Skip only pharokka | add `SKIP_PHAROKKA=1` |
| Skip pharokka+phold, start at FoldSeek (step 02) | add `SKIP_PHOLD=1` |
| Skip pharokka+phold+FoldSeek, start at compare (step 03) | add `SKIP_FOLDSEEK=1` |
| Skip through compare, start at curate (step 04) | add `SKIP_COMPARE=1` — the one to reach for when you've only edited curation rules and want to re-run `steps/04_curate.sh` (curate + build output, one job) against an existing run's FoldSeek + comparison output |

The three `SKIP_*` flags all need an **existing run** to resume from, so pair
them with `RESUME=1` or `PHAGEFACTOR_RUN_DIR=/…/run_00N`:

```bash
RESUME=1 SKIP_COMPARE=1 bash submit_all.sh
```

`submit_all.sh` enforces this: if `SKIP_PHOLD`/`SKIP_FOLDSEEK`/`SKIP_COMPARE`
is set without `RESUME=1` or `PHAGEFACTOR_RUN_DIR`, it refuses to start rather
than silently creating a new, empty run — sourcing your env file alone does
not set either of these, so it isn't enough by itself. It also checks that the
run you point it at actually has the output the flag expects to reuse (e.g.
`SKIP_COMPARE=1` needs that run's `03_comparison/comparison_per_gene.csv` to
already exist).

To resume from an **arbitrary** step in an existing run, submit that step and the
rest by hand with `--chdir` set to the run dir, e.g. from the merge onward:

```bash
source env_mymachine.sh
export PHAGEFACTOR_RUN_DIR=/…/run_002
R=$PHAGEFACTOR_RUN_DIR
J1c=$(sbatch  --parsable --chdir=$R steps/01c_merge.sh)
J02=$(sbatch  --parsable --chdir=$R --dependency=afterok:$J1c steps/02_foldseek_3di.sh)
# ... 03_compare.sh, 04_curate.sh (curates + builds output), then 05_phynteny.sh
```

## Step 05 (phynteny + synteny) on its own

Step 05 is optional and runs after step 04. `submit_all.sh` prints the exact
command; run it in a shell where you've `source`d your env file (so it picks up
the env + finds the run via `.current_run`):

```bash
source env_mymachine.sh
export PHAGEFACTOR_RUN_DIR=/…/run_002              # the run you want to extend
sbatch --chdir=$PHAGEFACTOR_RUN_DIR steps/05_phynteny.sh
```

## Installing on a cluster

Clone the repo, then create your `env_*.sh` on the cluster:

```bash
git clone https://github.com/TMs-code/phagefactor.git
```

A quick sanity check before submitting:

```bash
python -m py_compile scripts/*.py      # silent = OK
bash -n config.sh submit_all.sh                         # silent = OK
```
