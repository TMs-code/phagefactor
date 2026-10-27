# tests

`run_workflow_tests.sh` — end-to-end tests of the Snakemake workflow with mock tools
(no database, GPU or SLURM; ~2-3 min). See the header of the script and
[`docs/snakemake.md`](../docs/snakemake.md#5-tests).

`mock_bin/` — stand-ins for `pharokka.py`, `phold`, `foldseek`, `phynteny_transformer`
and `micromamba` that write the real tools' file layouts with deterministic synthetic
content. They exist only for these tests; never put them on PATH for a real run.
