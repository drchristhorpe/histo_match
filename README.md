# histo_match

Match an **MHC Class I allele sequence** to known allele numbers (e.g.
`HLA-A*02:01`), using [`histo_hmm`](https://github.com/drchristhorpe/histo_hmm)'s
profile HMMs to narrow the search down to a handful of candidate loci
before matching against ~14,200 bundled reference sequences.

It ships as:

- a Python library — `import histo_match`
- a CLI tool — `histo-match`
- a [Claude Code / Claude Desktop / Claude Science skill](skills/histo-match/SKILL.md)

Requires Python 3.14+.

## Install

```bash
uv sync                 # dev environment, from a checkout
uv tool install .       # install the `histo-match` CLI globally
# or
pip install .
```

## CLI

```bash
histo-match match --sequence "MAVMAPRTLVLLLSG..."
histo-match match --fasta sequences.fasta --top 10
histo-match match --fasta sequences.fasta --top 10 --json-file results.json
histo-match match --sequence "MAVMAPRTLVLLLSG..." --plain
```

`match` reports, per sequence:

- whether it's an **exact match** to a known reference sequence, or the
  best **approximate match** (Levenshtein identity) when it isn't
- the assigned allele: an `allele_slug` (e.g. `hla-a_02_01`), the full
  `gene_allele_name`/`protein_allele_name` (e.g. `HLA-A*02:01`), and the
  database `allele_id`
- the predicted locus and its HMM confidence, and a ranked table of the
  next-best candidate alleles

When multiple allele names share the identical reference sequence (or tie
on approximate identity), the one with the **lowest allele number** is
assigned — e.g. a sequence shared by `HLA-A*02:01` and `HLA-A*02:13`
resolves to `HLA-A*02:01`.

Options:

- `--sequence` — an amino-acid sequence. Repeat for multiple inputs.
- `--fasta` — match every sequence in a FASTA file.
- `--data-dir` — override the bundled reference data directory.
- `--candidate-loci` — number of HMM-ranked loci to search (default 3).
- `--min-locus-probability` — drop candidate loci below this HMM probability.
- `--top` — number of ranked matches to report when there's no exact match.
- `--plain` — compact machine-friendly stdout instead of Rich panels/tables.
- `--json-file PATH` — also write full results to a JSON file.

## Library

```python
from histo_match import AlleleMatcher

matcher = AlleleMatcher()
result = matcher.match("MAVMAPRTLVLLLSGALALTQTWA...")

print(result.exact_match)
print(result.predicted_locus, result.locus_confidence)
print(result.best_match.allele_slug)       # e.g. "hla-a_02_01"
print(result.best_match.gene_allele_name)  # e.g. "HLA-A*02:01"
print(result.matches)                      # ranked AlleleMatch list
```

`AlleleMatcher(data_dir=None, n_candidate_loci=3, min_locus_probability=0.0)`
loads the `histo_hmm` classifier once and reuses it — construct one
instance and call `.match()`/`.match_batch()` on it repeatedly, rather than
constructing a new one per sequence.

## Notes and limitations

- Reference sequences are matched after stripping whitespace and the `-`
  (unresolved residue) / `?` (unknown residue) placeholder characters used
  in the bundled data, uppercased — the same cleaning `histo_hmm` applies
  before scoring.
- Search-space reduction depends on the query's true locus appearing among
  the top `--candidate-loci` HMM predictions. For closely related paralogous
  loci (e.g. cattle `BoLA-1`/`BoLA-2`), the true locus doesn't always rank
  in the default top 3 — raise `--candidate-loci` for a broader (slower)
  search if exact matches are unexpectedly missed.
- Approximate matching ranks candidates by Levenshtein identity within the
  candidate loci only, not the full 251-locus reference set.

## Development

```bash
uv sync
uv run pytest
```

Tests use real, small loci already present in the bundled reference data
(e.g. `sla_6.json`) rather than synthetic fixtures.

See [docs/PLAN.md](docs/PLAN.md) for the design rationale and
[CHANGELOG.md](CHANGELOG.md) for release history.
