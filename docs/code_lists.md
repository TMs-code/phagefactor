# phageFACTor — where the curated word-lists / patterns live in the code

Companion to [`scoring_curation_DAG.svg`](scoring_curation_DAG.svg). Paths are
relative to the repo root.

Two files:

- **`scripts/lib/lexicon.py`** — all vocabulary (word-lists, patterns, the
  category↔keyword map). The definition site: `config.py`,
  `foldseek_scoring.py`, `utils.py`, `03_compare_annotations.py` and
  `04_curate_annotations.py` import from it and define no vocabulary of their own.
- **`scripts/lib/rules.py`** — the ordered rules that consume it (name upgrades,
  description → PHROG category, `short_name` generation).

`config.py` holds paths and thresholds; it re-exports three lexicon names
(`UNINFORMATIVE_STRINGS`, `GENERIC_WORDS`, `COMPLEMENTARY_CATEGORY_MAP`) so
`from config import ...` call sites keep working.

Search the symbol name rather than trusting a line number.

## lexicon.py, by section

| § | What it controls | Symbols |
|---|---|---|
| **A · informativeness** | the gate that drops most hits, plus the three ranking maluses | `UNINFORMATIVE_STRINGS` (exact blocklist) · `UNCHAR_PATTERN` (hypothetical / DUF / bare folds / accession-shaped names) · `GENERIC_DOMAIN` (×0.60, "family protein", "domain-containing", "repeat protein", bare "X-like", "motif") · `GENERIC_NAMES` (×0.60, UPF/Yabc/orf labels, TnpV) · `GENERIC_STRUCTURAL_PRODUCTS` (real but unspecific: `endolysin`, `major capsid protein`) · `GENERIC_WORDS` (tokens dropped before the Jaccard) · `GENERIC_DESCRIPTORS` (exact-membership set for the generic-vs-specific gate) |
| **B · taxonomy** | kingdom assignment without an NCBI taxonomy dependency | `BACTERIA_KEYWORDS` · `EUKARYOTE_KEYWORDS` · `ARCHAEA_KEYWORDS` · `VIRUS_KEYWORDS` · `EUKARYOTIC_DESC` (description-level euka signal) · `is_eukaryotic_taxname()` (word-boundary matcher shared by 03 and by 02d's parse-time filter) |
| **C · phage terms** | the boosts and the wrapper malus | `PHAGE_SPECIFIC` (×2.0) · `PHAGE_CONTEXT` (×1.5) · `GP_WRAPPER` (×0.40) · `PHAGE_WRAPPER_ONLY` (anchored: blocks a bare `gpNN` from winning) |
| **D · defense** | one vocabulary, two term groups | `DEFENSE_SYSTEM_TERMS` (defense in the broad sense: anti-phage / immune systems and counter-defense — restriction alleviation, anti-restriction, anti-CRISPR, immunity → category `defense`) · `DEFENSE_TA_TERMS` (toxin-antitoxin families → category `rules.TA_CATEGORY`) · `DEFENSE_ANY` = both = the FoldSeek defense flag (step 02) and the curation check `_desc_is_defense` (step 04) · `EUKA_RESCUE` (keeps Sir2/TIR/cGAS/Argonaute-type names that look eukaryotic). The defense/TA rows of `rules.CATEGORY_RULES` are built from the same two term groups. |
| **E · promiscuous folds** | fold ≠ function: flag + cap at BORDERLINE | `PROMISCUOUS_FOLDS` |
| **F · PDB titles** | "Crystal structure of X" → X; drop ORF/DARPin/SeMet titles | `PDB_CRYSTAL_PREFIX` · `PDB_TRAILING_QUAL` · `PDB_UNINFORMATIVE` · `PDB_TITLE` · `extract_pdb_description()` |
| **G · characterized domains** | the positive counterpart of E — domain-level folds that *do* carry a corroborated activity (SPOR, PIN, LysM, HNH, GIY-YIG…), each with activity + PMID | `RELEVANT_DOMAINS` · `relevant_domain()` |
| **H · curation gates** | the two deliberately different widths used in step 04 | `FS_GENERIC` (narrow) · `FS_GENERIC_BROAD` (wide, a superset of `PHAGE_WRAPPER_ONLY`) · `GENE_LABEL` (subset of `GENERIC_NAMES`; drives an upgrade, not a malus). The wrapper patterns are built from shared fragments, so they stay nested. |
| **I · functional equivalence** | category ↔ keyword map: the `complementary` agreement label, and (symmetrically) the synonym groups used by `_shares_function` | `COMPLEMENTARY_CATEGORY_MAP` |

Public API: `classify(desc) -> (tier, tags)` and `is_informative_desc(desc)`.
`tier` (0 uninformative / 1 low / 2 proper) matches the three levels of
[inPhormer](https://github.com/Mila-MP/inPhormer). It is read only by the
FoldSeek-only low-information gate and the top-3 majority vote; every other
predicate reads `tags`, because the callers need different widths and each width
drives a different consequence (drop vs low-information malus ×0.60 vs upgrade vs gate).
`classify` is cached per description, so repeated hits cost one evaluation.

## Ordered rules: `scripts/lib/rules.py`

The other half. A lexicon entry is a set member — order is irrelevant, and adding
one widens a test by exactly that one string. A rule is (pattern → action) and
**first match wins**, so a rule added at the top can shadow every rule below it.
Mixing the two in one file makes that difference invisible, so they are separate.

| Symbol | Used by | What it does |
|---|---|---|
| `UPGRADE_RULES` | step 04 | (desc pattern, top-3 pattern, replacement) — swap a generic winning name for a specific one already present in the top-3 hits. Supports `\1` backreferences. |
| `CATEGORY_RULES` · `category_from_text()` · `phrog_category()` · `PHROG_CANON` | steps 04, 05, 07 | description → PHROG category, one ordered map (specific phrases first, generic enzyme classes last). Step 04 uses it when PHold leaves the category unknown, step 05 for pharokka passthrough rows. |
| `make_short_name()` (`_HEAD_NOUNS` · `_DEFENSE_SHORT` · `_DECORATION` · `_LEADING` · `_GENE_SYMBOL`) | steps 05, 07 | `final_product` → a concise `short_name`. A pure post-pass: it cannot change *which* annotation was chosen. Manual exceptions go in `input/overrides.tsv`. |

## Elsewhere on purpose

| Symbol | File | Why it stays |
|---|---|---|
| `_RE_REG` · `_RE_INT` · `_RE_STR` · `_RE_REP` · `_RE_REG_EXCLUDE` | `scripts/lib/synteny_hint.py` | neighbourhood classes for the optional C1/Cro synteny hint; used nowhere else. |
| `_SUBDB_SOURCES` · `_SUBDB_FILE_RE` | `01p_merge_phold_proteins.py` · `03_compare_annotations.py` | phold sub-DB filenames, i.e. an interface to another tool. |

## Notes

- Section A is where *most* hits are removed (every hypothetical / DUF / unknown
  structural homolog).
- A dropped hit and a demoted hit are not the same thing: C's boosts and A's
  ×0.60 (low-information) / ×0.40 (wrapper) maluses are multiplicative on the ranking score only; E caps confidence;
  B and D set flags or demote one step.
- Eukaryotic hits are touched at three stages, from one list:
  `is_eukaryotic_taxname()` removes them *before* scoring (local search only),
  `EUKARYOTIC_DESC` / the same keyword list flag or demote them *after*, and
  `EUKA_RESCUE` brings back the ones that are meaningful in phage defence
  (Sir2, TIR, cGAS, Argonaute). Matching is on word boundaries, not substrings:
  *Campylobacter conciSUS* is not a pig.
