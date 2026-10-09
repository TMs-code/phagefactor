# Curation gate logic, phage-boost factor & naming rules

This is the scientific heart of phageFACTor (the **Annotate + Curate** in FACT).
It reconciles **Phold** (phage-DB structural/profile) and **custom FoldSeek**
(general structural) evidence per hypothetical gene. Implemented in
`03_compare_annotations.py` (compare) and `04_curate_annotations.py` (curate),
with scoring in `foldseek_scoring.py`.

**Every description word-list and regex lives in `lexicon.py`** — it is the single
definition site; the other modules import from it and define no vocabulary of their
own. Edit a term there once and every gate sees it. The lists are deliberately kept
*separate* rather than merged into one `is_generic()`, because each drives a
different consequence (drop / score malus / description upgrade / narrow gate);
`lexicon.py`'s docstring explains the widths. `classify(desc)` returns
`(tier, tags)` where tier is 0/1/2 to match inPhormer. The gates read `tags`;
`tier` is read only by the FoldSeek-only low-information gate (§6) and the top-3
majority vote (§3).

## 1. FoldSeek confidence tiers (`_compute_fs_confidence`)

| Tier | Condition (e = evalue, s = score) |
|---|---|
| **CONFIDENT** | `e ≤ 1e-3 AND s ≥ 200`, **or** top-3 agreement with (`e ≤ 1e-3` or `s ≥ 300`) |
| **GOOD** | `e ≤ 0.01` or `s ≥ 200`; or top-3 agreement with `e ≤ 0.05` |
| **BORDERLINE** | `e ≤ 0.1` or (`s ≥ 90` with `e ≤ 0.5`); or top-3 agreement at `e ≤ 0.1`/`s ≥ 90` |
| **WEAK** | anything below — excluded from auto-annotation |
| **NO_HIT** | no hit returned |

All thresholds are set in `config.yaml` (`foldseek.confident_*`, `good_*`, `borderline_*`).

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
| Generic suffix ("family protein", "domain-containing protein", "repeat protein", bare "X-like", "motif", …) and gene labels | **× 0.60** (`foldseek.low_info_malus`) |
| Wrapper ("gpNN", "phage protein", "Mu-like prophage …") | **× 0.40** |

Terms are matched as whole words ("holin" does not match "acetylcholine").

**Description upgrade (`_upgrade_description`).** The best hit's description may be
replaced by a name that occurs at least twice in the top-3 (majority vote), but only
when this does not override stronger evidence: a title, a eukaryote-specific name, a
promiscuous fold, a generic descriptor or a low-information label never replaces a
named best hit, and otherwise the majority's best raw score must be ≥ that of the
best hit. A generic, low-information or eukaryote-specific best hit is replaced by a
majority that names a function.

(The name also nods to the factor in phage**FACTor**.)

## 4. False-positive filters

- **Promiscuous folds** (`_PROMISCUOUS_FOLD_PATTERNS`): MBL/TLD/β-lactamase,
  glyoxalase II, tRNase Z, CPSF-73/100; eukaryotic motor proteins
  (kinesin/myosin/dynein); ERAD/ubiquitin ligase; nephrocystin; **Arc/HicB
  ribbon-helix-helix**. Flagged → routed to `needs_review` (Aravind 1999
  ; Daiyasu 2001). *Not* auto-rejected — real positives exist.
  The Arc/HicB ribbon-helix-helix case is demoted **only when BORDERLINE**; a
  CONFIDENT/GOOD hit is still reviewed.
- **PDB/structure-title detector** (`_is_pdb_title`): sentence-like pdb100 titles
  (descriptive structure-deposition strings rather than a function) are treated as
  uninformative.
- **Eukaryotic kingdom filter** (`03`, `_apply_eukaryotic_demote`): if best hit is
  Eukaryote and no Bacteria/Archaea/Virus appears in top-3, demote one tier
  (prevents Photosystem-I-style convergent-fold false positives). This same
  condition is also surfaced as a standalone boolean, `foldseek_euka_kingdom`
  (not just a tier change), specifically so `04_curate_annotations.py`'s
  `foldseek_only` case (no phold call at all -- see §6) can still catch a hit
  that got demoted down to BORDERLINE with a generic-sounding description that the
  description-keyword filter below would otherwise miss entirely.
- **Eukaryotic description filter** (`foldseek_euka_desc`): keyword fallback for
  DBs without embedded taxonomy (afdb-swissprot, pdb100) — matches organism
  names/eukaryote-specific gene families directly in the description text
  (`EUKARYOTIC_DESC` in `lexicon.py`). Weaker than the
  kingdom filter above: it only fires when the description text itself names an
  organism or a recognisably eukaryotic gene family, so a generic gene-symbol
  description with correct-but-unmentioned eukaryotic taxonomy needs the
  kingdom-based flag to be caught.

## 5. Agreement classification (step 03)

`strong` (Jaccard ≥ 0.35) · `partial` (≥ 0.08) · `complementary` (category matches
FS keyword via `COMPLEMENTARY_CATEGORY_MAP`, `lib/lexicon.py`) · `phold_only` · `foldseek_only` ·
`different` · `both_uninformative`.

## 6. Curation decision tree (step 04)

| Agreement | Action | Review? |
|---|---|---|
| strong / partial | auto-merge (prefer Phold unless FS is a strict superset) | no |
| complementary | auto-merge (FS if it is the direct functional match) | no |
| phold_only | use Phold | no |
| foldseek_only | use FoldSeek; **promiscuous/eukaryotic → review** | only if flagged |
| different | curation gate tree below, else flag | yes if unresolved |
| both_uninformative | "hypothetical protein" | no |

### Curation gate tree (grouped by function)

Sub-DB substitution runs upstream in `03_compare_annotations.py`, before
agreement is even classified — it is a standalone step, not part of the tree
below. Everything else is a pre-check at the top of the `different` branch in
`merge_annotations()` (`04_curate_annotations.py`): each is tried **in the
order it appears in the code**; the first one that matches auto-resolves the
gene, without ever reaching the relatedness-fallback gate at the end.

Gates are grouped by what KIND of decision they make, not by where they sit in
the file — a group is a documentation/lookup device (`grep '\[merge:'
04_curate_annotations.py` finds every merge gate), it does not reorder
anything. Only the **keep-phold** group has a real within-group priority
(broadest / most-confident phold-keeping check first); the other groups are
disjoint patterns where relative order doesn't matter.

**Standalone (upstream, outside the tree)**

| Gate | What it does |
|---|---|
| sub-DB substitution | (upstream, `03_compare_annotations.py`) copies a specialised sub-DB hit's structured name (VFDB/CARD/DefenseFinder/NetFlaX/ACR) into `phold_product` *before* agreement is classified, so downstream logic sees the real identity from the start. |

**Specialised sub-DB gate** — runs before the agreement cases, for genes whose PHold top hit is a NetFlaX / VFDB / CARD entry. These are reference collections, so the hit is weighed by its identity to the entry: *close* = `phold_fident` ≥ `curate.subdb_close_min_fident` (0.4; `phold_fident` = fraction of identical residues in PHold's alignment to the entry), *remote* = below. A FoldSeek call "names a function" when it is informative, not a title, not a narrow/broad generic phrase, not an exact generic descriptor (`Lipoprotein`, `Secreted protein`, …) and not a low-information name (see `fs-only` below).

| Situation | Result | `curation_action` |
|---|---|---|
| calls conflict (`different`), FoldSeek CONFIDENT + names a function, sub-DB remote | FoldSeek name `[keep-foldseek: subdb-yield]` | `auto_fs_over_weak_phold` |
| NetFlaX, FoldSeek GOOD/CONFIDENT + names a function (otherwise) | `"<FS name> (NetFlaX-predicted toxin\|antitoxin)"` `[merge: ta-functionalize]` | `auto_merge_TA_functionalised` |
| NetFlaX, FoldSeek generic / weak / absent | `"NetFlaX-predicted toxin\|antitoxin"` `[keep-phold: subdb-priority]` | `auto_phold_over_weak_fs` |
| VFDB / CARD otherwise | unchanged (Cases 1–5, `subdb-priority`) | unchanged |

NetFlaX naming: `NetFlaX-predicted toxin` / `NetFlaX-predicted antitoxin` when the NetFlaX record's `type` is toxin/antitoxin, otherwise `NetFlaX-predicted toxin-antitoxin system protein`. NetFlaX keeps the toxin-antitoxin category `rules.TA_CATEGORY` (*moron, auxiliary metabolic gene and host takeover*, PHROG convention) and is not forced to `defense`; the partner accession is not part of the name (it remains in `sub_db_tophits/`).

**`fs-only`** — checks on a FoldSeek-only call (Case 4):

| Gate | What it does |
|---|---|
| title-fallback | the description is a structure title (`… variant Delta+PHS L36R at cryogenic temperature`, `cA4-bound Cami1 in complex with 70S ribosome`) → the protein name inside the title (`Cami1`), else the first top-3 name that names a function, else no final annotation. Length alone (≥ 8 words) never demotes. |
| low-information | a FoldSeek-only name that is lexicon tier 1 (domain/family/repeat/motif suffix, gene label, e.g. `SHOCT domain-containing protein`, `XkdX family protein`, `YodL-like`) or a broad wrapper — with no characterised domain (`RELEVANT_DOMAINS`) and no function category inferable from it (`rules.CATEGORY_RULES`, so terminase / nuclease / DNA-binding labels keep their name) — is not used as the final annotation, at any confidence. Result: `hypothetical protein`, `both_uninformative`; the hit stays in `foldseek_description` / `foldseek_confidence` / `foldseek_top3` (and `fs_description` / `fs_top3` in the final table), the name in `curation_suggestion`. Generic-but-real descriptors (`Lipoprotein`, `DNA-binding protein`, `Secreted protein`) are not low-information and stay. |

**`keep-phold`** — phold's call wins outright, no merge, no review (priority order below):

| Gate | What it does |
|---|---|
| defense-priority | a confident (high/medium) Phold defense-system call is kept over a divergent FoldSeek fold; defense modules routinely reuse mobile-element folds, so disagreement here is expected and not a red flag. |
| atpase-family | Phold helicase/ABC-transporter/transporter + FoldSeek AAA(+)-ATPase/SMC are the same nucleotide-binding superfamily, not a disagreement. |
| weak-fs-guard | a BORDERLINE/WEAK or generic FoldSeek hit must not override or flag an informative Phold call; only CONFIDENT/GOOD FS competes. |
| transaldolase-keep | Phold "transaldolase" vs FS panB/hydroxymethyltransferase fold-alike → keep transaldolase. |
| suppress-fs-wrapper | suppress FS hits that are *only* a bare "(pro)phage protein" wrapper — an informative qualifier in front ("baseplate phage protein") is NOT suppressed. |
| top3-corroboration | scan *all* FS top-3 (hyphen-normalised) for corroboration; gated to `fs_conf == CONFIDENT`. Low-confidence Phold IS included — if Phold's call (even low conf) matches any FS top-3 entry, that cross-method agreement is honoured (e.g. Phold low "kinase" + FS top-3 "histidine kinase" → keep "kinase"). Per-hit FS confidence isn't in the top-3 string, so parity is enforced only at the gene level. |

**`keep-foldseek`** — FoldSeek's call wins outright, no merge, no review:

| Gate | What it does |
|---|---|
| confident-fs-relevant | a CONFIDENT/GOOD, non-promiscuous, **phage-relevant** FoldSeek call is taken over a weak Phold guess. "Phage-relevant" requires a named phage-specific/viral-context description, a same-host-genus hit (plausible moron/AMG), or a defense hit — a structurally confident but out-of-context FS hit does **not** qualify and falls through to review instead. |

**`merge`** — combine both calls into one annotation (disjoint patterns, order doesn't matter):

| Gate | What it does |
|---|---|
| ta-functionalize | Phold names a toxin-antitoxin role (NetFlaX/antitoxin/toxin/AbrB) and FoldSeek names the actual protein → merge as `"<FS name> (<role> protein)"`. |
| transposition-motor | Phold "DNA transposition"/transposase + FoldSeek AAA-ATPase/NTPase-KAP/P-loop → the FS hit is the transposition motor; merge `"DNA transposition (<motor>)"`. |
| morphogenesis | two structural/assembly terms (tail, head, portal, capsid, baseplate, …) on both sides → parts of the same morphogenesis module; merge `"<Phold> (<FoldSeek>)"`. |
| subdb-over-generic-fs | a specialised sub-DB hit (adopted by the sub-DB substitution step) beats a *generic* FoldSeek call even at low Phold confidence. |
| fs-specificity-upgrade | Phold gives a generic structural term and ≥2 FS top-3 hits converge on a specific named product → merge `"<generic> (<specific>)"`. |

**`specificity`** — prefer whichever side is genuinely specific (symmetric):

| Gate | What it does |
|---|---|
| generic-vs-specific | one side a pure generic descriptor, other names a specific protein/family → prefer the specific one (exact-membership check, never collides with `rules.UPGRADE_RULES`). |

**`fallback`** — nothing above matched:

| Gate | What it does |
|---|---|
| relatedness-fallback | nothing above resolved it. If Phold and FS (top-1 or any top-3) share a specific functional token, or Phold's category is complementary to FS, they aren't really divergent → `auto_merge_related`, no flag. Only genuinely unrelated disagreements reach `needs_review_divergent`. |

Each gate's return explanation is tagged `[group: slug]` (e.g. `[keep-phold:
defense-priority]`) so the source of an auto-curated call is traceable
directly from `curation_explanation` in the output table, and the
`curation_action` values are functional names.

### Characterized-domain post-process (`lexicon.py`, section G)
Runs *after* the full decision tree above, on every gene's final call
(`04_curate_annotations.py::_apply_relevant_domain`) — not gated to any
particular agreement branch or rule. It re-checks the resulting description
against a short list of structurally-recognisable, well-corroborated domains
(SPOR, PIN, LysM, HNH, GIY-YIG) via a domain-anchored regex ("`<token> domain
[-containing] protein`", not a specific protein name that merely contains the
token). A domain-level hit can never be more specific than a real functional
call, so: if Phold is informative, its product is kept with the domain
appended in brackets (`"lipoprotein (SPOR domain)"`); if Phold is
absent/uninformative, the domain's characterized activity is used instead
(`"ribonuclease (PIN domain)"`). Either outcome is always flagged for review
— domain-level evidence is relevant but not a confident final call.

### Source attribution: `both agree` vs `merged`
Both come from the Phold+FoldSeek path, but are distinguished for diagnostics:
- **`both agree`** — only when `agreement == strong` (high Jaccard, near-identical hit).
- **`merged`** — every other combined case (partial / complementary /
  suppress-fs-wrapper / fs-specificity-upgrade / generic-vs-specific /
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
  (the same curated map, reused symmetrically — *no new constant*; defined in
  `scripts/lib/lexicon.py` section I). New synonym
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

### Function-category inference (`rules.CATEGORY_RULES`)
When Phold leaves the category `unknown`, the description is matched against one
ordered keyword map (also used at output for pharokka passthrough rows). Specific
phrases come first (terminase/portal/capsid → *head and packaging*; tail, connector,
integration, lysis), then defense and toxin-antitoxin, transcription regulation,
DNA/RNA metabolism, host-derived metabolic enzymes (*moron, auxiliary metabolic gene
and host takeover*), and generic enzyme classes (*other*) last. Phold's own category
calls are never overwritten.

## 7. Specialised sub-DB integration

Phold often records only a generic placeholder when ACR/VFDB/CARD/NetFlaX/
DefenseFinder wins. `03` reads `sub_db_tophits/*.tsv` directly to recover the
structured name (gene symbol, system, ARO name) and corrects the category to
`defense` for DefenseFinder/ACR. NetFlaX is named by role only
(`NetFlaX-predicted toxin` / `antitoxin`), keeps Phold's TA category and is
resolved by the specialised sub-DB gate in `04` (section on gates above).

## 8. Naming rules → `short_name` (`lib/rules.py`, `make_short_name`)

Deterministic post-pass (does not change which annotation was chosen):
1. Sub-DB `Name (Category): desc` → gene symbol (RecN, HlyA, …) when symbol-like.
2. Slash alternatives → first (primary) call.
3. Strip decoration: `(Fragment)`, `(EC …)`, "from bacteriophage X", "domain-containing protein" → "protein", leading Phage/Putative.
4. Trailing gene symbol after ≥3 descriptor words → the symbol.
5. Long enzyme names → head-noun class (aldolase, kinase, …).
6. Else trim to ≤5 words.

## 9. Manual overrides (`input/overrides.tsv`, applied in `05_build_output.py`)

`input/overrides.tsv` (`match_type, match, final_product, final_function,
short_name, annotation_source, note`) is the **escape hatch** for irreducibly
per-gene decisions that should not become regex rules: forcing the FoldSeek call
(`choose-FS`), specific merges, gene-symbol `short_name`s, etc. It applies in
`05_build_output.py` (run as part of `steps/04_curate.sh`) rather than only in
`07_integrate.py` — so a hand-curated call lands in the deliverable table **and
clears the review flag** for that gene. The `annotation_source` column lets an
override declare the evidence category it should carry in the final table.

## 10. Phynteny + synteny integration (`07_integrate.py`)

Runs after `06_phynteny.py` (Phynteny Transformer) — both scripts are invoked by the
single optional `steps/05_phynteny.sh` SLURM step. Merges two synteny-based signals
into the final table. Phynteny IDs are **positional** (e.g. `<record>_0` = the first
CDS in that GBK record), so they are mapped to `locus_tag` by CDS order within each
record (exact; works in genome and protein mode).

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
