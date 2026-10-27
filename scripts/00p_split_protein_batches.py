#!/usr/bin/env python3
# =============================================================================
# 00p_split_protein_batches.py — split protein-mode input into fixed-size
#                                 phold batches
#
# A fixed batch size (default 50 proteins, workflow setting protein_batch_size)
# keeps batch boundaries stable: editing one protein changes one batch, so only
# that batch re-runs phold. The last batch may be partial.
#
# Input  : <fasta dir>/*.{faa,fa,fasta}  (one file per prophage/source;
#          PROTEIN_FASTA_DIR, default <input>/fasta)
# Output : <run>/split/all_proteins_combined.faa      (">PROPH:ID" headers)
#          <run>/split/protein_batches/batch_NNNN.faa
#          <run>/split/protein_batch_list.txt         (one batch name per line)
# =============================================================================

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from config import INPUT_DIR, PROTEIN_COMBINED_FAA as COMBINED_FAA, \
    PROTEIN_BATCH_DIR as BATCH_DIR, PROTEIN_BATCH_LIST as BATCH_LIST

# Set by the Snakefile (protein_batch_size; pangenome mode points PROTEIN_FASTA_DIR
# at the translated FASTA).
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
    # whatsoever - not phage-like". Catch it here, before phold runs, so a
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
        log("    - input_mode: auto   (auto-detect translates it: 'pangenome'), or")
        log("    - --config input_mode=pangenome   (explicit), or")
        log("    - translate yourself: scripts/00n_translate_pangenome.py")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 2: split into fixed-size batches (last batch may be partial)
    # ------------------------------------------------------------------
    n_batches_expected = (n_total + BATCH_SIZE - 1) // BATCH_SIZE

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


if __name__ == "__main__":
    main()
