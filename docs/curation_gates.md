# Curation gate logic, phage-boost factor & naming rules

This is the scientific heart of phageFACTor (the **Assign + Curate** in FACT).
It reconciles **Phold** (phage-DB structural/profile) and **custom FoldSeek**
(general structural) evidence per hypothetical gene. Implemented in
`03_compare_annotations.py` (compare) and `04_curate_annotations.py` (curate),
with scoring in `foldseek_scoring.py`.

## 1. FoldSeek confidence tiers (`_compute_fs_confidence`)

| Tier | Condition (e = evalue, s = score) |
|---|---|
| **CONFIDENT** | `e ≤ 1e-3 AND s ≥ 200`, **or** top-3 agreement with (`e ≤ 1e-3` or `s ≥ 300`) |
| **GOOD** | `e ≤ 0.01` or `s ≥ 200`; or top-3 agreement with `e ≤ 0.05` |
| **BORDERLINE** | `e ≤ 0.1` or (`s ≥ 90` with `e ≤ 0.5`); or top-3 agreement at `e ≤ 0.1`/`s ≥ 90` |
| **WEAK** | anything below — excluded from auto-annotation |
| **NO_HIT** | no hit returned |

Score-only CONFIDENT was deliberately removed: high structural-alignment scores
are common for promiscuous folds (see §4).

## 2. Quality gate (step 02)

A hit is kept if `evalue < FOLDSEEK_EVALUE_MAX (0.1)` **OR**
`score ≥ FOLDSEEK_SCORE_OVERRIDE (200)` (rescue). Per-residue ProstT5 3Di tokens
below `PROSTT5_MASK_THRESHOLD (25/100)` are masked to `*` before search.

## 3. Phage-boost factor (`_phage_boost_factor`)

Composite ranking score = `foldseek_score × boost`:

| Match in description | Boost |
|---|---|
| Named phage protein (terminase, capsid, portal, gpNN, holin, integrase, …) | **× 2.00** |
| Generic phage/viral context | **× 1.50** |
| Generic suffix ("family protein", "domain-containing protein", …) | **× 0.75** |

(The name also nods to the factor in phage**FACTor**.)

## 4. False-positive filters

- **Promiscuous folds** (`_PROMISCUOUS_FOLD_PATTERNS`): MBL/TLD/β-lactamase,
  glyoxalase II, tRNase Z, CPSF-73/100; eukaryotic motor proteins
  (kinesin/myosin/dynein); ERAD/ubiquitin ligase; nephrocystin; **Arc/HicB
  ribbon-helix-helix**. Flagged → routed to `needs_review` (Aravind 1999 PMID
  11471255; Daiyasu 2001 PMID 11513844). *Not* auto-rejected — real positives exist.
  The Arc/HicB ribbon-helix-helix case is demoted **only when BORDERLINE**; a
  CONFIDENT/GOOD hit is still reviewed.
- **PDB/structure-title detector** (`_is_pdb_title`): sentence-like pdb100 titles
  (descriptive structure-deposition strings rather than a function) are treated as
  uninformative.
- **Eukaryotic kingdom filter** (`03`, `_apply_eukaryotic_demote`): if best hit is
  Eukaryote and no Bacteria/Archaea/Virus appears in top-3, demote one tier
  (prevents Photosystem-I-style convergent-fold false positives).
- **Eukaryotic description filter**: keyword fallback for DBs without embedded
  taxonomy (afdb-swissprot, pdb100).

## 5. Agreement classification (step 03)

`strong` (Jaccard ≥ 0.35) · `partial` (≥ 0.08) · `complementary` (category matches
FS keyword via `COMPLEMENTARY_CATEGORY_MAP`) · `phold_only` · `foldseek_only` ·
`different` · `both_uninformative`.

## 6. Curation decision tree (step 04)

| Agreement | Action | Review? |
|---|---|---|
| strong / partial | auto-merge (prefer Phold unless FS is a strict superset) | no |
| complementary | auto-merge (FS if it is the direct functional match) | no |
| phold_only | use Phold | no |
| foldseek_only | use FoldSeek; **promiscuous/eukaryotic → review** | only if flagged |
| different | Rules A–E pre-checks, else flag | yes if unresolved |
| both_uninformative | "hypothetical protein" | no |

### Rules A–E (pre-checks on the `different` branch, order B → A → C → E)
- **Rule 0** (sub-DB low-conf): a specialised sub-DB hit (VFDB/CARD/T6SS…) beats a
  *generic* FoldSeek call even at low Phold confidence.
- **Rule B**: suppress FS hits that are *only* a "(pro)phage protein" wrapper — but
  keep informative qualifiers ("baseplate phage protein" is NOT suppressed).
- **Rule A**: Phold gives a generic structural term and ≥2 FS top-3 hits converge on
  a specific named product → merge `"<generic> (<specific>)"`.
- **Rule C**: scan *all* FS top-3 (hyphen-normalised) for corroboration; gated to
  `fs_conf == CONFIDENT`. **Low-confidence Phold IS included** — if Phold's call (even
  low conf) matches any FS top-3 entry, that cross-method agreement is honoured
  (e.g. Phold low "kinase" + FS top-3 "histidine kinase" → keep "kinase"). Per-hit FS
  confidence is not in the top-3 string, so parity is enforced only at the gene level.
- **Rule E**: one side a pure generic descriptor, other names a specific
  protein/family → prefer specific (exact-membership, never collides with `_UPGRADE_RULES`).

### Source attribution: `both agree` vs `merged`
Both come from the Phold+FoldSeek path, but are distinguished for diagnostics:
- **`both agree`** — only when `agreement == strong` (high Jaccard, near-identical hit).
- **`merged`** — every other combined case (partial / complementary / Rule-B/A/E /
  different-resolved): the two sources were *combined*, not identical.

Plus `phold`, `foldseek`, `no_hit`. Surfaced in the final table as the evidence category
(`structural-only` = FoldSeek-only).

### Relatedness check (`_shares_function`)
Decides whether the Phold and FoldSeek calls describe the *same* function. Beyond the
Jaccard/keyword tests, it now also treats two names as related via:
- **(a) split-token substring + shared prefix ≥ 6.** Catches compound or typo'd enzyme
  names (e.g. `metallo-protease` ~ `metallopeptidase`; a misspelled
  `deoxyribosyltransferase`) that token-set overlap alone would miss.
- **(b) synonym groups** sourced from the **`COMPLEMENTARY_CATEGORY_MAP` value-lists**
  (the same curated map, reused symmetrically — *no new constant*). New synonym
  members in that map: lysis `+ glycosidase / glucosaminidase / phosphodiester`;
  head & packaging `+ protease / peptidase / metallopeptidase / metalloprotease`;
  transcription `+ fur / ferric / regulation`; DNA-metabolism `+ soj`; moron
  `+ oxidoreductase / redox / fe-s / fes / ferredoxin`.
- **(c) ATPase superfamily** rule also covers **ABC-transporter / SMC**.

When Phold and FoldSeek are found "related", the single clean **Phold** name is kept
instead of emitting a verbose `"phold / foldseek"` string.

### Defense gate
A defense flag forces `needs_review` **only when the *chosen* Phold/FS call is
itself a defense annotation**.

### Function-category inference (`_infer_function_cat`, `_FUNC_CAT_RULES`)
When Phold leaves the category `unknown`, the description is matched against keyword
rules to assign a PHROG category. Notably, host-derived **auxiliary metabolic / host-
takeover (moron)** genes are recognised by enzyme/role keywords: `ribosom*`,
`zinc-uptake`/`zur`/`fur-family`, `kinase`, `oxidoreductase`, `dehydrogenase`,
`transaldolase`, `aldolase`, `racemase`, `reductase`, `permease`, `synthase`,
`isomerase`, `epimerase`, `mutase`, `phosphatase`, `sialyltransferase` → *moron,
auxiliary metabolic gene and host takeover*. Phold's own category calls are never
overwritten (the rule only fires on `unknown`).

## 7. Specialised sub-DB integration

Phold often records only a generic placeholder when ACR/VFDB/CARD/NetFlaX/
DefenseFinder wins. `03` reads `sub_db_tophits/*.tsv` directly to recover the
structured name (gene symbol, system, ARO name) and corrects the category to
`defense` for DefenseFinder/NetFlaX/ACR.

## 8. Naming rules → `short_name` (`naming_rules.py`)

Deterministic post-pass (does not change which annotation was chosen):
1. Sub-DB `Name (Category): desc` → gene symbol (RecN, HlyA, …) when symbol-like.
2. Slash alternatives → first (primary) call.
3. Strip decoration: `(Fragment)`, `(EC …)`, "from bacteriophage X", "domain-containing protein" → "protein", leading Phage/Putative.
4. Trailing gene symbol after ≥3 descriptor words → the symbol.
5. Long enzyme names → head-noun class (aldolase, kinase, …).
6. Else trim to ≤5 words.

## 9. Manual overrides (`input/overrides.tsv`, applied in step 05)

`input/overrides.tsv` (`match_type, match, final_product, final_function,
short_name, annotation_source, note`) is the **escape hatch** for irreducibly
per-gene decisions that should not become regex rules: forcing the FoldSeek call
(`choose-FS`), specific merges, gene-symbol `short_name`s, etc. It now applies in
**step 05** (build output) rather than only step 07 — so a hand-curated call lands in
the deliverable table **and clears the review flag** for that gene. The new
`annotation_source` column lets an override declare the evidence category it should
carry in the final table.

## 10. Phynteny + synteny integration (step 07, `07_integrate.py`)

Runs after step 06 (Phynteny Transformer). Merges two synteny-based signals into the
final table. Phynteny IDs are **positional** (`phiNP_0` = Nth CDS in the GBK record),
so they are mapped to `locus_tag` by CDS order within each record (exact; works in
genome and protein mode).

- **Phynteny → `final_function` (priority).** When Phynteny is confident
  (`≥ threshold`, default 0.8) and names a real PHROG category, it **overwrites**
  `final_function` (most agree; some fill gaps left `unknown`; the few that differ →
  Phynteny wins, e.g. Zur → "transcription regulation"). Pipeline category is kept only
  where Phynteny is unknown/other or below threshold. `phynteny_probability` is kept as
  a far-right column; no separate phynteny category column (avoid column sprawl).
- **C1/Cro → `short_name`.** `synteny_hint.py` finds the divergent immunity switch
  (adjacent opposite-strand regulators, gap < 1.5 kb) and assigns **CI/C1 repressor**
  (faces the lysogeny/integrase side, larger) vs **Cro** (faces the structural/lytic
  side, small) by transcription direction + size. A `[likely]` call (orientation AND
  size agree) **replaces** `short_name`; tentative calls are left as a hint only.
- **`synteny_hint` column** also tags integration boundaries, lysis cassettes and the
  late/structural operon. It is indicative (replaces nothing else) and is written into
  the final GBK as CDS `/note` qualifiers.
