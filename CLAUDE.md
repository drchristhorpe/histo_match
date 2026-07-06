# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`histo_match` matches an MHC Class I protein sequence to known allele
numbers (e.g. `HLA-A*02:01`). It uses
[`histo_hmm`](https://github.com/drchristhorpe/histo_hmm) to predict the
sequence's most likely locus/loci, then matches against the reference
sequences bundled for those candidate loci only — this is the search-space
reduction the tool exists for. It's a Python library, a Click-based CLI
(`histo-match`), and a Claude skill wrapping that CLI
(`skills/histo-match/`). See [README.md](README.md) for user-facing usage
and [docs/PLAN.md](docs/PLAN.md) for the design rationale.

## Environment

- Python 3.14, managed with `uv`. Use `uv sync`, `uv run <cmd>`, `uv run pytest`.
- Don't invoke a bare `python`/`pip` — always go through `uv run` /
  `uv add` so the lockfile stays authoritative.
- `histo_hmm` is a git dependency (`[tool.uv.sources]` in `pyproject.toml`
  points at `github.com/drchristhorpe/histo_hmm`, branch `master`), not a
  local path — it doesn't assume a sibling checkout exists.

## Layout

```
src/histo_match/
  core.py            # AlleleMatch, MatchResult, AlleleMatcher (the matching algorithm)
  reference_data.py  # per-locus JSON loading/caching, allele_sort_key, slugify_allele_name, clean_sequence
  cli.py             # Click CLI (entry point: histo-match)
  data/cytoplasmic_sequences/*.json   # bundled reference data, ships with the package
tests/
  test_reference_data.py
  test_core.py
  test_cli.py
skills/histo-match/SKILL.md
```

## Key invariants — don't break these

- `AlleleMatcher.__init__` loads the `histo_hmm` classifier **once**;
  `.match()`/`.match_batch()` must reuse it, not reconstruct per call.
- The reference data's own `canonical_allele` field is unreliable — 77 of
  14,180 bundled entries have it as an empty dict despite a populated
  `alleles` list. **Never read `canonical_allele` from the raw JSON.**
  `LocusReferenceData` always self-derives the canonical allele as the
  lowest-numbered entry in `alleles` (`allele_sort_key` on
  `gene_allele_name`) — this is also the user's explicit rule for
  resolving ties, not just a data-cleanliness workaround.
- Query and reference sequences must both be cleaned with
  `reference_data.clean_sequence()` (uppercase, strip whitespace, strip
  `-`/`?`) before any comparison. 1,557 of the 14,180 bundled sequences
  contain a literal `-`, marking a genuinely unresolved residue at that
  position (not an MSA gap) — comparing an uncleaned reference key against
  a cleaned query silently breaks exact-match lookup.
- Exact matching is an O(1) dict lookup per candidate locus, tried in HMM-rank
  order — stop at the first hit, don't fall through to approximate matching
  once one is found. Approximate matching only runs when no candidate locus
  has an exact hit.
- Approximate-match ties, and the ordering used to build
  `LocusReferenceData.ordered_sequences`, must resolve to the **lowest
  allele number** — this is why `ordered_items()` iterates ascending by
  `allele_sort_key` and why the final sort in `AlleleMatcher.match()` is a
  stable sort keyed only on `-identity` (ties keep their ascending-order
  position).
- Search-space reduction is not guaranteed to include the true locus —
  closely related paralogous loci (e.g. `bola_1`/`bola_2`) can outrank each
  other in `histo_hmm`'s `top_loci`. This is expected, not a bug to "fix"
  by hardcoding locus groups; `n_candidate_loci`/`min_locus_probability`
  are the user-facing knobs.

## Testing

- `uv run pytest` — no network access needed; the bundled reference data
  under `src/histo_match/data/` is used directly as test fixtures (e.g.
  `sla_6.json`, chosen because `histo_hmm` classifies its sequences as
  `sla_6` at rank 1 with full confidence, so tests aren't flaky against the
  HMM's own ranking behavior).
- If you add a fixture sequence for an exact-match test, verify with the
  real classifier first which locus it actually ranks under
  `n_candidate_loci` — don't assume a sequence's source locus file name is
  what `histo_hmm` will predict for it.
- `AlleleMatcher()` construction (~0.5s, loads 251 profile HMMs) is done
  once per test module via a `scope="module"` fixture in `test_core.py` —
  don't construct a fresh one per test function.

## Scope

The CLI has one command, `match`. Don't add further subcommands or output
formats without checking with the user first.
