#!/usr/bin/env python3
# =============================================================================
# 00p_split_protein_batches.py — split protein-mode input into fixed-size
#                                 phold batches for SLURM array submission
#
# Replaces the old "concatenate everything, run phold as one big GPU job"
# approach in 01p_phold_proteins.sh with fixed 50-protein batches, run as a
# CPU-default SLURM array job (01p_phold_proteins_array.sh) — mirrors the
# genome-mode 01_phold_array.sh pattern.
#
# Why fixed-size batches (not "split into N chunks"):
#   User decision (2026-06): "lets have 50-prot batches for step01 (dont
#   divide by 4, i want it to be reproducible for next datasets, so 50 is a
#   good number, and its ok if the last batch is not full)". A fixed batch
#   size means batch boundaries — and therefore array-task counts, job names,
#   and `.done` markers — are stable and reproducible across datasets of any
#   size, whereas "divide into N chunks" would shift every boundary whenever
#   the input count changes.
#
# This script is meant to run INLINE (plain `python3 ...`, not via sbatch)
# from submit_all.sh, *before* the array job is submitted — the array size
# (`--array=0-$((N_BATCHES-1))`) can only be computed once batches exist.
# It is fast (pure I/O, no model inference) so this is safe.
#
# Input  : PROTEIN_FASTA_DIR/*.{faa,fa,fasta}  (one file per prophage/source)
# Output : INPUT_DIR/all_proteins_combined.faa      (prefixed concatenation)
#          INPUT_DIR/protein_batches/batch_NNNN.faa (fixed PROTEIN_BATCH_SIZE seqs)
#          INPUT_DIR/protein_batch_list.txt         (one batch name per line,
#                                                     drives the array job)
#
# Re-run safety: if protein_batch_list.txt already exists and its batch count
# matches the number of batch files on disk, splitting is skipped (idempotent
# — safe to call from submit_all.sh on every invocation).
# =============================================================================

import os
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))

# INPUT_DIR comes from the python config module (config.py -> config.py,
# the same passthrough every other scripts/*.py uses) so the path always
# matches what the rest of the pipeline resolves — no risk of drifting from a
# bash-only value.
try:
    from config import INPUT_DIR as _CFG_INPUT_DIR
    INPUT_DIR = Path(_CFG_INPUT_DIR)
except ImportError:
    INPUT_DIR = Path(os.environ.get("INPUT_DIR", str(PROJECT_DIR / "input")))

# PROTEIN_BATCH_SIZE and PROTEIN_FASTA_DIR are bash-only config.sh values
# (config.py has no equivalent) — submit_all.sh exports them before
# calling this script. Defaults here mirror config.sh (PROTEIN_BATCH_SIZE=50,
# PROTEIN_FASTA_DIR=${INPUT_DIR}/fasta) so standalone runs behave identically.
FASTA_DIR = Path(os.environ.get("PROTEIN_FASTA_DIR", str(INPUT_DIR / "fasta")))
BATCH_SIZE = int(os.environ.get("PROTEIN_BATCH_SIZE", "50"))

COMBINED_FAA = INPUT_DIR / "all_proteins_combined.faa"
BATCH_DIR = INPUT_DIR / "protein_batches"
BATCH_LIST = INPUT_DIR / "protein_batch_list.txt"


def log(msg):
    print(f"[00p_split] {msg}", flush=True)


def _iter_fasta_records(path):
    """Yield (header_line_without_gt, [sequence_lines]) for a FASTA file."""
    header, seq_lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header, seq_lines = line[1:], []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, seq_lines


def main():
    if not FASTA_DIR.is_dir():
        log(f"ERROR: PROTEIN_FASTA_DIR not found: {FASTA_DIR}")
        sys.exit(1)

    protein_files = sorted(
        p for p in FASTA_DIR.iterdir()
        if p.suffix.lower() in (".faa", ".fa", ".fasta") and p.is_file()
    )
    if not protein_files:
        log(f"ERROR: no protein FASTA files (.faa/.fa/.fasta) found in {FASTA_DIR}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 1: concatenate with per-file prophage prefix (">ID" -> ">PROPH:ID")
    # Mirrors the `sed "s/^>/>$_PROPH:/"` logic previously inline in
    # 01p_phold_proteins.sh, so downstream gene-ID parsing (split on ":")
    # remains unchanged.
    # ------------------------------------------------------------------
    log(f"Concatenating {len(protein_files)} FASTA file(s) -> {COMBINED_FAA}")
    records = []  # list of (prefixed_header, seq_lines)
    for f in protein_files:
        proph = f.stem
        if proph.endswith("_proteins"):
            proph = proph[: -len("_proteins")]
        n_seqs = 0
        for header, seq_lines in _iter_fasta_records(f):
            records.append((f"{proph}:{header}", seq_lines))
            n_seqs += 1
        log(f"  {proph}: {n_seqs} sequences from {f.name}")

    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(COMBINED_FAA, "w") as out:
        for header, seq_lines in records:
            out.write(f">{header}\n")
            for sl in seq_lines:
                out.write(sl + "\n")

    n_total = len(records)
    log(f"Total protein sequences: {n_total}")

    # ------------------------------------------------------------------
    # Step 2: split into fixed-size batches (last batch may be partial)
    # ------------------------------------------------------------------
    n_batches_expected = (n_total + BATCH_SIZE - 1) // BATCH_SIZE

    if BATCH_LIST.is_file():
        existing = [l.strip() for l in BATCH_LIST.read_text().splitlines() if l.strip()]
        existing_files = sorted(BATCH_DIR.glob("batch_*.faa")) if BATCH_DIR.is_dir() else []
        if len(existing) == n_batches_expected == len(existing_files):
            log(f"protein_batch_list.txt already up to date "
                f"({len(existing)} batches of <= {BATCH_SIZE} proteins) — skipping split")
            log(f"Batch list: {BATCH_LIST}")
            return

    BATCH_DIR.mkdir(parents=True, exist_ok=True)
    # Clear stale batch files from a previous run with a different batch count
    for stale in BATCH_DIR.glob("batch_*.faa"):
        stale.unlink()

    batch_names = []
    for i in range(0, n_total, BATCH_SIZE):
        batch_idx = i // BATCH_SIZE
        batch_name = f"batch_{batch_idx:04d}"
        batch_path = BATCH_DIR / f"{batch_name}.faa"
        chunk = records[i:i + BATCH_SIZE]
        with open(batch_path, "w") as out:
            for header, seq_lines in chunk:
                out.write(f">{header}\n")
                for sl in seq_lines:
                    out.write(sl + "\n")
        batch_names.append(batch_name)
        log(f"  {batch_name}: {len(chunk)} proteins -> {batch_path}")

    BATCH_LIST.write_text("\n".join(batch_names) + "\n")
    log(f"Wrote {len(batch_names)} batches (fixed size {BATCH_SIZE}, last batch "
        f"has {n_total - (len(batch_names) - 1) * BATCH_SIZE} proteins)")
    log(f"Batch list -> {BATCH_LIST}")
    log(f"Array job should be submitted as: --array=0-{len(batch_names) - 1}")


if __name__ == "__main__":
    main()
