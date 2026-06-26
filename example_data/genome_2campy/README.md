# Example 1 — genome mode (2 Campylobacter prophages)

The simplest mode: input is just the nucleotide genome FASTA of each prophage.
Pharokka (Phanotate) calls the CDS, then phold + custom FoldSeek annotate the
hypotheticals, and the FACT curation layer assigns final functions.

```
genome_2campy/input/
├── fasta/
│   ├── CJIE1.fasta        # one nucleotide FASTA per prophage (header = >NAME)
│   └── CJIE2.fasta
└── prophage_list.txt      # the prophage basenames, one per line
```

Two representative *Campylobacter jejuni* prophages (~33 kb each), CJIE1 and CJIE2.

## Run it
```bash
# from the repo root, point input/ at this example
cp -r example_data/genome_2campy/input/* input/
bash submit_all.sh                 # auto-detects "genome" from input/fasta/
#   ...or the no-database path:
SEARCH_MODE=webapi bash submit_all.sh
```

Expected: `04_output/final_annotations_table.csv` with one row per CDS, the
hypotheticals resolved by phold/FoldSeek, and `04_output/curation/review_suggested.csv`
listing the few genes flagged for a human look. See `../README.md` for the column
reference and `../run_smoke_test.sh` for an automated check.
