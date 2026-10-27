#!/usr/bin/env python3
"""Mock pharokka.py: naive forward-strand ORF calling -> <prefix>.gbk + CDS TSV."""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _mocklib import rng
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

ap = argparse.ArgumentParser()
ap.add_argument("-i"); ap.add_argument("-o"); ap.add_argument("-d"); ap.add_argument("-t")
ap.add_argument("-p"); ap.add_argument("-f", action="store_true")
a = ap.parse_args()
out = Path(a.o); out.mkdir(parents=True, exist_ok=True)
rec = next(SeqIO.parse(a.i, "fasta"))
seq = str(rec.seq).upper()
feats, rows, i, n = [], [], 0, 0
stops = {"TAA", "TAG", "TGA"}
while i < len(seq) - 3:
    if seq[i:i + 3] == "ATG":
        j = i
        while j < len(seq) - 3 and seq[j:j + 3] not in stops:
            j += 3
        if j - i >= 300 and j < len(seq) - 3:
            n += 1
            lt = f"{a.p}_CDS_{n:04d}"
            prot = str(Seq(seq[i:j]).translate())
            annot = "hypothetical protein" if rng("ph", lt).random() < 0.7 else "tail protein"
            feats.append(SeqFeature(FeatureLocation(i, j + 3, strand=1), type="CDS",
                                    qualifiers={"locus_tag": [lt], "product": [annot],
                                                "translation": [prot]}))
            rows.append((lt, i + 1, j + 3, "+", annot))
            i = j + 3
            continue
    i += 1
gb = SeqRecord(Seq(seq), id=a.p, name=a.p, description=f"{a.p} mock pharokka",
               annotations={"molecule_type": "DNA"}, features=feats)
SeqIO.write([gb], out / f"{a.p}.gbk", "genbank")
with open(out / f"{a.p}_cds_final_merged_output.tsv", "w") as fh:
    fh.write("gene\tstart\tstop\tframe\tcontig\tannot\n")
    for lt, s, e, st, an in rows:
        fh.write(f"{lt}\t{s}\t{e}\t{st}\t{a.p}\t{an}\n")
print(f"mock pharokka: {n} CDS -> {out}")
