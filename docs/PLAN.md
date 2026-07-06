# histo_match design plan

## 1. Purpose

Given a query MHC Class I protein sequence, identify which known allele(s) it
corresponds to (e.g. resolve a sequence to `HLA-A*02:01`). A brute-force
search against all ~14,200 reference sequences across 251 loci is wasteful.
The sibling project [`histo_hmm`](https://github.com/drchristhorpe/histo_hmm)
trains one profile HMM per locus and can classify a query sequence's most
likely locus/loci in one call — this is the search-space reduction
`histo_match` relies on: use `histo_hmm` to shortlist candidate loci, then
only search reference sequences within those loci.

## 2. Reference data

`data/cytoplasmic_sequences/*.json` (251 files, one per locus, ~14,200
sequences total) is the same data `histo_hmm` was trained on. Each file maps
a reference amino-acid sequence to its allele metadata:

```json
{
  "<sequence>": {
    "alleles": [
      {"gene_allele_name": "...", "id": "...", "locus": "...",
       "protein_allele_name": "...", "source": "..."},
      ...
    ],
    "canonical_allele": { ... }
  }
}
```

Locus/file-stem names match `histo_hmm`'s model/locus names exactly (its
`manifest.json` classes list equals this repo's filename stems, e.g.
`hla_a`, `mamu_b11l`), so a predicted locus maps directly to `<locus>.json`.

Where a sequence is shared by multiple allele names, `canonical_allele` is
already the lowest-numbered allele among them (verified by inspection of
`hla_a.json`, e.g. `HLA-A*01:01:05` chosen over `:15`/`:26`/`:51`/...).
However, 77 of the 14,180 bundled entries (across several loci) have an
*empty* `canonical_allele` (`{}`) despite a populated `alleles` list — so
`histo_match` never trusts the source field directly. It always
self-derives the canonical allele as the lowest-numbered entry in
`alleles` (natural/numeric sort on `gene_allele_name`), which happens to
agree with `canonical_allele` everywhere it's populated, and is correct
by construction everywhere it isn't.

The data is packaged inside the installed library at
`src/histo_match/data/cytoplasmic_sequences/` (mirrors how `histo_hmm`
bundles its trained models under `src/histo_hmm/models/`), not left at the
repo root, so it ships with the wheel.

## 3. Dependency on histo_hmm

`histo_hmm` isn't published to PyPI, so it's added as a `uv` git source
pointing directly at GitHub (portable — doesn't assume a sibling checkout
exists on every machine):

```toml
dependencies = [..., "histo-hmm"]
[tool.uv.sources]
histo-hmm = { git = "https://github.com/drchristhorpe/histo_hmm.git", branch = "master" }
```

`histo_hmm.MHCClassIClassifier` is instantiated once (loading all 251
profile HMMs is the expensive part) and reused across calls. Its
`.classify(sequence)` returns `is_class_i`, `top_loci` (locus, probability
pairs), `best_score`, and `region_start`/`region_end` (for sequences
embedded in larger constructs).

## 4. Matching algorithm

1. Clean the query sequence (strip whitespace, uppercase, strip `-`/`?`
   placeholder characters — the same convention `histo_hmm` uses before
   training/scoring). Reference sequence keys are cleaned the same way at
   load time: 1,557 of the 14,180 bundled reference sequences (across 53
   loci) contain a literal `-` marking an unresolved/missing residue at
   that position (not an alignment gap), so exact-match lookup must compare
   like-for-like cleaned strings on both sides, not raw query vs. raw
   dict key.
2. Classify via `histo_hmm` to get `is_class_i`, `top_loci`,
   `region_start`/`region_end`.
3. Trim the sequence to `[region_start:region_end]`.
4. For each candidate locus in `top_loci` (bounded by `n_candidate_loci` /
   `min_locus_probability`), lazily load and cache that locus's reference
   JSON:
   - Try an O(1) dict lookup of the trimmed sequence. On a hit: this is an
     exact match. Take the entry's `canonical_allele` as the best match and
     its full `alleles` list as all tied matches. Stop — no need to check
     lower-ranked candidate loci (a dict key match is definitionally the
     best possible result).
   - If no candidate locus has an exact hit, fall back to approximate
     matching using [rapidfuzz's `Levenshtein`](https://github.com/rapidfuzz/Levenshtein)
     (`ratio`/`distance`, C-optimized). Within each candidate locus, iterate
     reference sequences ordered by ascending allele number (falling back
     to plain string comparison for non-numeric nomenclature, e.g.
     `Rano-E*u`) and score each against the query. Track the best (highest
     ratio) result across all candidate loci, keeping the first-seen
     (lowest allele number) entry on exact ties. Report the top-N as
     ranked matches.
5. If nothing scores usefully (e.g. `is_class_i` is false and no candidate
   locus yields a reasonable match), return an empty result rather than a
   forced guess.

**Ordering rationale:** processing candidates in ascending allele-number
order and stopping at the first hit is a deliberate design choice, not just
an optimization. Exact matching is already O(1) via dict lookup, so
early-exit's real payoff is at the *locus* level (skip lower-ranked
candidate loci once an exact hit is found) and as the natural tie-break
order for approximate matches (first-seen-wins on a distance tie, which is
always the lowest allele number given ascending iteration).

**`allele_slug`:** the match output's primary identifier is a derived
`allele_slug`, not the database `id` field. Lowercase the allele name and
replace `*`/`:` with `_`, derived from `gene_allele_name`, e.g.
`HLA-A*02:01` → `hla-a_02_01`, `Rano-E*u` → `rano-e_u`. The full
`gene_allele_name`/`protein_allele_name`/`allele_id`/`source` are retained
for completeness, but the slug is what CLI output and summaries lead with.

## 5. Library API

```python
from histo_match import AlleleMatcher

matcher = AlleleMatcher()
result = matcher.match("MAVMAPRTLVLLLSGALALTQTWA...")

print(result.exact_match)
print(result.best_match.allele_slug)
print(result.matches)
```

- `AlleleMatch`: `allele_slug`, `gene_allele_name`, `protein_allele_name`,
  `allele_id`, `locus`, `source`, `identity`, `matched_sequence`.
- `MatchResult`: `query_sequence`, `cleaned_sequence`, `is_class_i`,
  `predicted_locus`, `locus_confidence`, `candidate_loci`, `region_start`,
  `region_end`, `exact_match`, `best_match`, `matches`.
- `AlleleMatcher(data_dir=None, n_candidate_loci=3, min_locus_probability=0.0)`
  — loads the `histo_hmm` classifier once; `.match(sequence, n_top=5)` and
  `.match_batch(sequences, n_top=5)`.

## 6. CLI

```
histo-match match --sequence "MAVMAP..." [--top N] [--json-file out.json] [--plain]
histo-match match --fasta sequences.fasta [--top N] [--json-file out.json] [--plain]
```

Mirrors `histo-hmm classify`'s flag shape for consistency between the two
tools. FASTA parsing via Biopython's `SeqIO`. Human-readable output by
default; `--json-file` writes structured results; `--plain` for
machine-friendly stdout.

## 7. Claude skill

`skills/histo-match/SKILL.md`, same frontmatter shape as `histo_com`'s and
`histo_hmm`'s skills (`name` + trigger-oriented `description`). Documents
`histo-match match` usage and cross-references the `histo-hmm` skill for
locus-only classification — no code coupling, just documentation.

## 8. Package layout

```
histo_match/
  .gitignore
  .python-version
  pyproject.toml
  README.md
  CLAUDE.md
  CHANGELOG.md
  docs/
    PLAN.md
  src/histo_match/
    __init__.py
    core.py               # AlleleMatcher, MatchResult, AlleleMatch
    reference_data.py     # lazy per-locus JSON loading/caching, sort key, slugify()
    cli.py                # Click CLI entry point `histo-match`
    py.typed
    data/cytoplasmic_sequences/*.json
  skills/histo-match/SKILL.md
  tests/
    test_reference_data.py
    test_core.py
    test_cli.py
  tmp/                    # gitignored scratch dir for manual smoke-test runs
```

## 9. Testing plan

Real, small loci already bundled in the data (e.g. `bola_1.json`, 17
sequences) are used directly as fixtures — fast and deterministic, no need
for synthetic data. Cover: exact match including the multi-allele tie →
lowest number via `canonical_allele`; approximate match via `Levenshtein`
ranking with the ascending-allele-number tie-break; locus narrowing through
the real `histo_hmm` classifier; and CLI invocation via
`click.testing.CliRunner`.

## 10. Workflow

1. Write this plan, commit it.
2. Scaffold project with `uv init`/`uv add`; move `data/` into
   `src/histo_match/data/`.
3. Implement `reference_data.py`, `core.py`, `cli.py`.
4. Write and run tests.
5. Write README.md, CLAUDE.md, CHANGELOG.md, and the Claude skill.
6. Manually exercise exact match, approximate (single-mismatch) match, and
   a construct-embedded/partial sequence end-to-end via the CLI.
7. Commit, on approval.
