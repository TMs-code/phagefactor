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
#   A fixed batch size (default 50 proteins, config.sh PROTEIN_BATCH_SIZE)
#   means batch boundaries — and therefore array-task counts, job names, and
#   `.done` markers — are stable and reproducible across datasets of any
#   size, with the last batch simply left partial. "Divide into N chunks"
#   would instead shift every boundary whenever the input count changes.
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

PROJECT_DIR = next((d for d in [Path(__file__).resolve().parent,
                                *Path(__file__).resolve().parents]
                    if (d / "config" / "config.yaml").is_file()),
                   Path(__file__).resolve().parent.parent)
sys.path.insert(0, str(PROJECT_DIR))
# scripts/lib/ holds the shared modules (config, utils, ...); repo-root config/
# only holds config.yaml.
sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

# INPUT_DIR comes from the python config module (config.py -> config.py,
# the same passthrough every other scripts/*.py uses) so the path always
# matches what the rest of the pipeline resolves — no risk of drifting from a
# bash-only value.
# INPUT_DIR (data in) + the generated-artifact paths (under the run dir's split/)
# both come from config.py so they always match the rest of the pipeline.
try:
    from config import (INPUT_DIR as _CFG_INPUT_DIR, SPLIT_DIR as _CFG_SPLIT_DIR,
                        PROTEIN_COMBINED_FAA as _CFG_COMBINED,
                        PROTEIN_BATCH_DIR as _CFG_BATCH_DIR,
                        PROTEIN_BATCH_LIST as _CFG_BATCH_LIST)
    INPUT_DIR = Path(_CFG_INPUT_DIR)
    COMBINED_FAA = Path(_CFG_COMBINED)
    BATCH_DIR = Path(_CFG_BATCH_DIR)
    BATCH_LIST = Path(_CFG_BATCH_LIST)
except ImportError:
    INPUT_DIR = Path(os.environ.get("PHAGEFACTOR_INPUT", str(PROJECT_DIR / "input")))
    _SPLIT = Path(os.environ.get("SPLIT_DIR", str(PROJECT_DIR / "runs" / "run" / "split")))
    COMBINED_FAA = _SPLIT / "all_proteins_combined.faa"
    BATCH_DIR = _SPLIT / "protein_batches"
    BATCH_LIST = _SPLIT / "protein_batch_list.txt"

# PROTEIN_BATCH_SIZE and PROTEIN_FASTA_DIR are bash-only config.sh values
# (config.py has no equivalent) — submit_all.sh exports them before
# calling this script. Defaults here mirror config.sh (PROTEIN_BATCH_SIZE=50,
# PROTEIN_FASTA_DIR=${INPUT_DIR}/fasta) so standalone runs behave identically.
FASTA_DIR = Path(os.environ.get("PROTEIN_FASTA_DIR", str(INPUT_DIR / "fasta")))
BATCH_SIZE = int(os.environ.get("PROTEIN_BATCH_SIZE", "50"))


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
    # Downstream gene-ID parsing splits on ":".
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

    COMBINED_FAA.parent.mkdir(parents=True, exist_ok=True)
    with open(COMBINED_FAA, "w") as out:
        for header, seq_lines in records:
            out.write(f">{header}\n")
            for sl in seq_lines:
                out.write(sl + "\n")

    n_total = len(records)
    log(f"Total protein sequences: {n_total}")

    # ------------------------------------------------------------------
    # Sanity guard: refuse NUCLEOTIDE input in protein mode.
    # phold proteins-mode does NOT translate; if these are DNA (ATCGN only)
    # ProstT5 still emits (meaningless) 3Di and FoldSeek then finds "no hits
    # whatsoever - not phage-like". Catch it here, before any sbatch, so a
    # whole phold run is not wasted.
    # ------------------------------------------------------------------
    sample = records[: min(50, len(records))]
    tot = naa = 0
    for _h, _sl in sample:
        s = "".join(_sl).strip().upper()
        tot += len(s)
        naa += sum(1 for c in s if c in "EFHIKLMPQRSVWY")  # protein-only residues
    if tot and (naa / tot) < 0.005:
        log(f"ERROR: fasta/ looks like NUCLEOTIDE, not protein "
            f"({naa}/{tot} protein-specific residues across a {len(sample)}-seq sample).")
        log("  phold protein mode does NOT translate DNA. Do one of:")
        log("    - unset INPUT_MODE   (let auto-detect translate it: 'pangenome'), or")
        log("    - INPUT_MODE=pangenome bash submit_all.sh   (explicit), or")
        log("    - translate yourself: scripts/00n_translate_pangenome.py")
        sys.exit(1)

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
