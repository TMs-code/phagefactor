#!/usr/bin/env python3
"""
detect_input_type.py
====================
Auto-detect whether the FASTAs in input/fasta/ are:
  - nucleotide phage genomes  →  "genome"   (Mode 1: run Pharokka → Phold → FoldSeek)
  - protein sequences          →  "protein"  (Mode P: skip Pharokka, run Phold proteins mode)

Detection rules (applied in order):
  1. Length heuristic: if ANY sequence is > LENGTH_THRESHOLD nucleotides
     → genome  (even the smallest known phages are ~5 kb)
  2. Amino-acid heuristic: if ALL sequences contain at least one protein-specific
     amino acid (E/D/F/H/I/K/L/M/P/Q/R/S/V/W/Y) → protein
  3. Tie-breaker: if sequences are short but look like nucleotide (only ATCGN/X),
     and no protein-specific amino acids are present → default to genome.

The detected mode is printed to stdout AND written to --output (default
input/.input_type) so submit_all.sh can cache the result across runs.

Usage:
  python scripts/lib/detect_input_type.py                        # auto-finds input/fasta/
  python scripts/lib/detect_input_type.py --fasta-dir path/     # explicit FASTA dir
  python scripts/lib/detect_input_type.py --redetect            # ignore cached .input_type
  python scripts/lib/detect_input_type.py --force protein       # override

Exits 0 with "genome" or "protein" on stdout.
Exits 1 on error (no FASTA files found, ambiguous content).
"""

import argparse
import sys
import os
from pathlib import Path

LENGTH_THRESHOLD = 5000          # sequences longer than this → genome
MIN_PROTEIN_FRACTION = 0.02      # at least 2% of chars must be protein-specific AA

# --- bulk-nucleotide-CDS (pan-genome) guard ---------------------------------
# Genome mode expects a HANDFUL of genome-length (>5 kb) nucleotide FASTAs, one
# per prophage. A clustered pan-genome / pan-proteome arrives instead as MANY
# short, in-frame nucleotide CDS in a single file. Rule 3 below would silently
# call that "genome" and Phanotate-recall every tiny record — wrong. When we see
# many short nucleotide-only records we STOP with guidance instead of guessing.
PANGENOME_MIN_RECORDS = int(os.environ.get("PANGENOME_MIN_RECORDS", "20"))
PANGENOME_MAX_MEAN_LEN = int(os.environ.get("PANGENOME_MAX_MEAN_LEN", "3000"))

# Amino acids found exclusively in proteins (absent in nucleotide sequences)
PROTEIN_ONLY_AA = set("EFHIKLMPQRSVWYeefhiklmpqrsvwy")


def read_fasta_sequences(fasta_path: Path) -> list[tuple[str, str]]:
    """
    Parse a FASTA file and return list of (header, sequence) tuples.
    Only reads the first 10,000 characters of each sequence for speed.
    """
    seqs = []
    current_header = None
    current_seq = []
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    seqs.append((current_header, "".join(current_seq)))
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                if len("".join(current_seq)) < 10000:
                    current_seq.append(line.strip())
    if current_header is not None:
        seqs.append((current_header, "".join(current_seq)))
    return seqs


def classify_sequence(seq: str) -> str:
    """
    Classify a single sequence as 'genome', 'protein', or 'ambiguous'.
    """
    s = seq.upper().replace("-", "").replace("*", "")
    if not s:
        return "ambiguous"

    # Rule 1: length
    if len(s) > LENGTH_THRESHOLD:
        return "genome"

    # Rule 2: protein-specific amino acids
    n_protein_aa = sum(1 for c in s if c in PROTEIN_ONLY_AA)
    if n_protein_aa / len(s) >= MIN_PROTEIN_FRACTION:
        return "protein"

    # Rule 3: all nucleotide chars → probably genome (short contig)
    nucleotide_chars = set("ATCGNatcgnXx")
    if all(c in nucleotide_chars for c in s):
        return "genome"

    return "ambiguous"


def detect(fasta_dir: Path, verbose: bool = False) -> str:
    """
    Analyse all FASTAs in fasta_dir and return 'genome' or 'protein'.
    Raises SystemExit(1) if the content is genuinely ambiguous or no files found.
    """
    fasta_files = sorted(
        list(fasta_dir.glob("*.fasta")) +
        list(fasta_dir.glob("*.fa")) +
        list(fasta_dir.glob("*.faa")) +
        list(fasta_dir.glob("*.fna"))
    )
    if not fasta_files:
        print(f"ERROR: No FASTA files found in {fasta_dir}", file=sys.stderr)
        sys.exit(1)

    votes = {"genome": 0, "protein": 0, "ambiguous": 0}
    nuc_short_lens = []   # lengths of short, nucleotide-only records (pan-genome tell)
    for fpath in fasta_files:
        seqs = read_fasta_sequences(fpath)
        for header, seq in seqs:
            verdict = classify_sequence(seq)
            votes[verdict] += 1
            s = seq.upper().replace("-", "").replace("*", "")
            if s and len(s) <= LENGTH_THRESHOLD and all(c in "ATCGNX" for c in s):
                nuc_short_lens.append(len(s))
            if verbose:
                print(f"  {fpath.name}::{header}: len={len(seq)} → {verdict}",
                      file=sys.stderr)

    if verbose:
        print(f"  Votes: {votes}", file=sys.stderr)

    # --- bulk-nucleotide-CDS (pan-genome) detection -------------------------
    # MANY short nucleotide-only records = a pan-genome / bulk CDS set, NOT a set
    # of prophage genomes. Return the dedicated "pangenome" mode (instead of
    # silently picking genome mode and Phanotate-recalling each record);
    # submit_all.sh then translates it to protein and runs protein mode.
    if (votes["protein"] == 0
            and len(nuc_short_lens) >= PANGENOME_MIN_RECORDS
            and (sum(nuc_short_lens) / len(nuc_short_lens)) <= PANGENOME_MAX_MEAN_LEN):
        mean_len = sum(nuc_short_lens) / len(nuc_short_lens)
        print(
            f"INFO: input/fasta/ is a BULK NUCLEOTIDE CDS set / pan-genome "
            f"({len(nuc_short_lens)} short nucleotide records, mean {mean_len:.0f} nt), "
            "not per-prophage genomes -> 'pangenome' mode. submit_all.sh translates it "
            "to protein (00n_translate_pangenome.py) then runs protein mode. "
            "Override with --force genome|protein or INPUT_MODE=.",
            file=sys.stderr,
        )
        return "pangenome"

    # Genome takes priority: any genome-sized sequence → genome mode
    if votes["genome"] > 0:
        return "genome"

    # All protein votes → protein mode
    if votes["protein"] > 0 and votes["genome"] == 0:
        return "protein"

    # All ambiguous (very short sequences with only ATCGN) → genome by default
    if votes["ambiguous"] > 0 and votes["genome"] == 0 and votes["protein"] == 0:
        if verbose:
            print(
                "  WARNING: All sequences ambiguous (short, ATCGN-only). "
                "Defaulting to genome mode.",
                file=sys.stderr,
            )
        return "genome"

    print("ERROR: Mixed genome/protein content in input/fasta/. "
          "Keep genomes and proteins in separate runs.", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Detect whether input/fasta/ contains phage genomes or protein sequences."
    )
    parser.add_argument(
        "--fasta-dir", type=Path, default=None,
        help="Path to FASTA directory (default: input/fasta/ relative to script)"
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="File to write result (default: <fasta_dir>/../.input_type)"
    )
    parser.add_argument(
        "--redetect", action="store_true",
        help="Re-run detection even if .input_type cache file exists"
    )
    parser.add_argument(
        "--force", choices=["genome", "protein"],
        help="Override detection and force a specific mode"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print per-sequence classification details to stderr"
    )
    args = parser.parse_args()

    # Resolve FASTA dir
    if args.fasta_dir:
        fasta_dir = args.fasta_dir.resolve()
    else:
        # Default: input/fasta/ relative to the phageFACTor root
        script_dir = Path(__file__).parent
        fasta_dir = (script_dir.parent / "input" / "fasta").resolve()

    if not fasta_dir.exists():
        print(f"ERROR: FASTA directory not found: {fasta_dir}", file=sys.stderr)
        sys.exit(1)

    # Resolve output file
    output_file = args.output if args.output else (fasta_dir.parent / ".input_type")

    # Check cache
    if not args.redetect and not args.force and output_file.exists():
        cached = output_file.read_text().strip()
        if cached in ("genome", "protein", "pangenome"):
            print(cached)
            return

    # Force override
    if args.force:
        mode = args.force
    else:
        mode = detect(fasta_dir, verbose=args.verbose)

    # Write cache
    output_file.write_text(mode + "\n")

    print(mode)


if __name__ == "__main__":
    main()
