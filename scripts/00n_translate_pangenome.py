#!/usr/bin/env python3
# =============================================================================
# 00n_translate_pangenome.py — prepare a bulk NUCLEOTIDE CDS set (e.g. a Roary
#                              pan-genome) for phageFACTor PROTEIN mode.
#
# WHY THIS STEP EXISTS
# --------------------
# phageFACTor protein mode (00p_split -> 01p phold-proteins) is *protein-in*: it
# never translates. A pan-genome / clustered pan-proteome delivered as
# nucleotide CDS (Roary representative-sequence FASTA, one gene per record,
# in-frame, starting near a start codon) therefore needs a one-shot translation
# BEFORE 00p_split runs. This is the "pangenome" sub-mode of protein mode
# (docs/modes.md, README "Example 3").
#
# WHAT IT DOES
# ------------
#   1. Reads a nucleotide CDS FASTA (each record = one gene, in frame).
#   2. Translates with a configurable codon table (default 11 = bacterial).
#      - trailing stop codon is stripped; a mid-sequence stop is a WARNING
#        (record is skipped by default — likely a frameshift / wrong table).
#   3. RENAMES every record to a clean, filesystem/tool-safe id
#      (<prefix>NNNNN, default prefix "pang_") so the dotted/pipe-laden source
#      ids (e.g. "CCIE1_AAAGWQ01.gff.1", "fig|123.4.peg.5") cannot be mangled
#      by phold/foldseek id handling.
#   4. Writes:
#        <out_faa>                      the translated protein FASTA (protein mode input)
#        <id_map>  (CSV)                clean_id, orig_id, group, nt_len, aa_len, note
#      The id map is the join key to recover the biological group name
#      (Roary "group_XXXX", the 2nd header token) after annotation — see
#      08_pangenome_reference.py.
#
# STANDALONE USAGE (no cluster / no config needed)
#   python scripts/00n_translate_pangenome.py \
#       --in  pangenome_nt.fa \
#       --out-faa input/fasta/pangenome_proteins.faa \
#       --id-map  pangenome_id_map.csv \
#       --table 11 --prefix pang_
#
# It writes input/.input_type = "protein" next to the fasta/ dir automatically
# (unless --no-input-type), so submit_all.sh's auto-detector is bypassed.
# =============================================================================

import argparse
import sys
from pathlib import Path


def log(msg):
    print(f"[00n_translate] {msg}", flush=True)


def iter_fasta(path):
    """Yield (id, description, sequence) for each record in a FASTA file."""
    hdr, desc, seq = None, "", []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, desc, "".join(seq)
                parts = line[1:].split(None, 1)
                hdr = parts[0] if parts else ""
                desc = parts[1] if len(parts) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
    if hdr is not None:
        yield hdr, desc, "".join(seq)


def translate(nt, table):
    """Translate a nucleotide string. Returns (protein, note).
    Uses Biopython if available, else a built-in table-11/standard codon map."""
    nt = nt.upper().replace("U", "T")
    trimmed = ""
    if len(nt) % 3 != 0:
        trimmed = f"trimmed_{len(nt) % 3}nt"
        nt = nt[: len(nt) - (len(nt) % 3)]
    try:
        from Bio.Seq import Seq
        prot = str(Seq(nt).translate(table=table))
    except Exception:
        prot = _translate_fallback(nt, table)
    note_bits = [b for b in (trimmed,) if b]
    # strip a single trailing stop
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot, ";".join(note_bits)


# Minimal fallback (bacterial table 11 == standard code for translation;
# alt start codons are translated to their actual residue, not forced to M).
_CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def _translate_fallback(nt, table):
    if table not in (11, 1):
        log(f"WARNING: fallback translator only supports tables 1/11; got {table}. "
            "Install biopython for other tables. Using standard code.")
    out = []
    for i in range(0, len(nt), 3):
        out.append(_CODON.get(nt[i:i + 3], "X"))
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description="Translate a bulk nucleotide CDS "
                                 "pan-genome into protein-mode input for phageFACTor.")
    ap.add_argument("--in", dest="infile", required=True, type=Path,
                    help="input nucleotide CDS FASTA (one gene per record, in frame)")
    ap.add_argument("--out-faa", required=True, type=Path,
                    help="output protein FASTA (place under input/fasta/)")
    ap.add_argument("--id-map", required=True, type=Path,
                    help="output CSV mapping clean_id -> orig_id/group/lengths")
    ap.add_argument("--table", type=int, default=11,
                    help="NCBI codon table (default 11 = bacterial)")
    ap.add_argument("--prefix", default="pang_",
                    help="clean-id prefix (default 'pang_')")
    ap.add_argument("--group-token", type=int, default=1,
                    help="0-based header token holding the biological group name "
                         "(Roary puts 'group_XXXX' at token 1; use -1 to disable)")
    ap.add_argument("--keep-internal-stop", action="store_true",
                    help="keep (don't skip) records with an internal stop codon")
    ap.add_argument("--no-input-type", action="store_true",
                    help="do not write <fasta>/../.input_type = protein")
    args = ap.parse_args()

    if not args.infile.exists():
        log(f"ERROR: input not found: {args.infile}")
        sys.exit(1)

    args.out_faa.parent.mkdir(parents=True, exist_ok=True)
    args.id_map.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_out = n_skipped = 0
    width = 5
    rows = []
    with open(args.out_faa, "w") as out:
        for hdr, desc, nt in iter_fasta(args.infile):
            n_in += 1
            group = ""
            if args.group_token >= 0:
                toks = ((hdr + " " + desc).split())
                if len(toks) > args.group_token:
                    group = toks[args.group_token]
            prot, note = translate(nt, args.table)
            if "*" in prot:
                if args.keep_internal_stop:
                    note = (note + ";" if note else "") + "internal_stop_kept"
                    prot = prot.replace("*", "X")
                else:
                    log(f"WARNING: internal stop in {hdr} — skipping "
                        "(likely frameshift or wrong codon table; "
                        "use --keep-internal-stop to force)")
                    n_skipped += 1
                    continue
            clean_id = f"{args.prefix}{n_out + 1:0{width}d}"
            out.write(f">{clean_id}\n{prot}\n")
            rows.append((clean_id, hdr, group, str(len(nt)), str(len(prot)), note))
            n_out += 1

    with open(args.id_map, "w", encoding="utf-8") as mh:
        mh.write("clean_id,orig_id,group,nt_len,aa_len,note\n")
        for r in rows:
            # quote orig_id/group defensively (may contain commas/pipes)
            cid, oid, grp, ntl, aal, note = r
            mh.write(f'{cid},"{oid}","{grp}",{ntl},{aal},"{note}"\n')

    if not args.no_input_type:
        itype = args.out_faa.parent.parent / ".input_type"
        itype.write_text("protein\n")
        log(f"Wrote {itype} = protein (bypasses auto-detection)")

    log(f"Read {n_in} nucleotide records; wrote {n_out} proteins "
        f"({n_skipped} skipped).")
    log(f"  protein FASTA -> {args.out_faa}")
    log(f"  id map        -> {args.id_map}")
    if n_out == 0:
        log("ERROR: no proteins written — check codon table / input frame.")
        sys.exit(1)


if __name__ == "__main__":
    main()
