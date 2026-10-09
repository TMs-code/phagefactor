# Architecture

Two diagrams. Open the SVGs directly — they are plain text and render in GitHub.

## 1. Pipeline — [`pipeline_DAG.svg`](pipeline_DAG.svg)

Three input modes (genome, protein, pangenome) converging on one shared spine:

```
00c pharokka → 01 phold → 02 FoldSeek → 03 compare → 04 curate → 05 output
                                          ↑                    → 06 phynteny → 07 integrate
```

Step 02 forks: `02d` searches local FoldSeek databases, `02w` uses the public
`search.foldseek.com` API. Both produce identically shaped `best_hit.csv` /
`top3.csv`, so everything downstream is unaffected by which one ran.

**`02w` is experimental.** The web API accepts *structures* only, so it needs an
ESMFold folding step first (`scripts/lib/esmfold_cache.py`), and it is subject to
server-side rate limits and timeouts. Use `02d` for anything you intend to
publish.

**Next:** integrating phageFACTor's curation (steps 03–04) into
[FoldFirstAskLater](https://github.com/hannelorelongin/FoldFirstAskLater) (FFAL).

## 2. Scoring and curation — [`scoring_curation_DAG.svg`](scoring_curation_DAG.svg)

Where a FoldSeek hit is filtered, scored, ranked and finally named. The
superscript numbers refer to the word-lists catalogued in
[`code_lists.md`](code_lists.md). Prose version of the same logic:
[`curation_gates.md`](curation_gates.md).

Two things the diagram is meant to make obvious:

1. **The informativeness gate removes most hits**, before any scoring happens.
2. **Trusted calls are never re-annotated.** Only pharokka-hypothetical genes
   enter curation at all; a sequence-homology call is passed through untouched.