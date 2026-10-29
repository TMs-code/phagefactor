# Configuration

phageFACTor reads two files from the `config/` folder of the directory it is
launched from (the *project dir*):

| File | Holds | Edit before the first run |
|---|---|---|
| `config.yaml` | the science: database locations, FoldSeek tiers, curation thresholds, host genus | `databases.db_root`, `project.host_genus` |
| `workflow.yaml` | what to run and where: input/search mode, Phynteny on/off, SLURM partition/QOS, resources per step | optional (every key defaults to a sensible value) |

Full reference: [`docs/configure.md`](../docs/configure.md).

## Input

Put the sequences in `input/fasta/` next to `config/`:

- **genome mode:** one nucleotide FASTA per prophage, plus `input/prophage_list.txt`
  with one prophage name (FASTA file stem) per line;
- **protein mode:** one or more protein FASTA files (`.faa`);
- **pangenome mode:** one nucleotide FASTA of gene sequences.

The mode is detected from the FASTA files (`input_mode: auto`). Details:
[`docs/modes.md`](../docs/modes.md).

## Databases

The databases (Pharokka, Phold, FoldSeek: AFDB-SwissProt, PDB100, optionally
AFDB50) are not downloaded by the workflow. Download them once per machine and set
`databases.db_root` in `config/config.yaml` to their folder:

```bash
conda create -n phagefactor -c conda-forge -c bioconda phagefactor   # provides the helper
conda run -n phagefactor phagefactor databases /path/to/db_dir --all
```

[`docs/databases.md`](../docs/databases.md) lists the sizes and how to download
only a subset (the FoldSeek search can also run on the FoldSeek web server:
`search_mode: webapi`).

## Running as a Snakemake module

Deployed with [snakedeploy](https://snakedeploy.readthedocs.io), the workflow
runs with plain Snakemake:

```bash
snakedeploy deploy-workflow https://github.com/TMs-code/phagefactor . --tag v1.2.0
snakemake --cores all --sdm conda                     # this machine
snakemake --sdm conda --executor slurm --jobs 50      # SLURM
```

The `config` handed to the module may hold both files' keys: the sections of
`config.yaml` (`project`, `paths`, `databases`, `foldseek`, `compare`, `curate`,
`phynteny`, `precomputed`) are merged over `config/config.yaml`, every other key
is a `workflow.yaml` setting. `--config phynteny=false` switches Phynteny off.

With `--sdm conda`, Snakemake builds `envs/phagefactor.yaml` and
`envs/phynteny.yaml` itself; the first Phynteny job downloads its models into
the env (needs internet). The Phold database must then be in `<db_root>/phold_db`
(or `databases.phold_db`): a Phold database installed inside another env is not
seen. Without `--sdm conda`, jobs run in the named
environments set by `tools_env` / `phynteny_env` (see
[`docs/snakemake.md`](../docs/snakemake.md)).
