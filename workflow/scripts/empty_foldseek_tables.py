#!/usr/bin/env python3
"""search_mode=none: write well-formed EMPTY best_hit.csv / top3.csv.

Step 03 requires both files; with no FoldSeek evidence the curation still runs on
pharokka + phold alone. Same columns as the real tables.
Usage: empty_foldseek_tables.py <best_hit.csv> <top3.csv>
"""
import sys
from pathlib import Path

import pandas as pd

best, top3 = (Path(a) for a in sys.argv[1:3])
best.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(columns=["gene", "accession", "description", "taxname", "pident", "evalue",
                      "score", "composite_score", "qcov_aa", "same_host", "defense_flag",
                      "informative_hit_found", "aa_length", "foldseek_subdb"]
             ).to_csv(best, index=False)
pd.DataFrame(columns=["gene", "rank", "description", "accession", "evalue", "score"]
             ).to_csv(top3, index=False)
print(f"wrote empty {best} and {top3} (search_mode=none)")
