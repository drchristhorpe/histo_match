# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-07-06

### Added

- `histo_match` Python library: `AlleleMatcher` class that uses
  `histo_hmm`'s profile HMMs to narrow a query sequence's candidate loci,
  then matches it against bundled reference sequences — an exact dict
  lookup first, falling back to Levenshtein-ranked approximate matching.
- `AlleleMatch`/`MatchResult` dataclasses, including a derived
  `allele_slug` (e.g. `hla-a_02_01`) as the primary output identifier.
- Ties (multiple allele names sharing one reference sequence, or equal
  approximate-match identity) resolve to the **lowest allele number**,
  self-derived per locus rather than trusting the source data's
  `canonical_allele` field, which is empty for 77 of 14,180 bundled entries.
- `histo-match` CLI (Click-based) with a `match` command:
  `--sequence`/`--fasta` input, `--candidate-loci`,
  `--min-locus-probability`, `--top`, `--plain`, `--json-file`.
- Bundled reference data (`src/histo_match/data/cytoplasmic_sequences/`,
  251 loci, ~14,200 sequences) ships inside the installed package.
- `histo_hmm` added as a git dependency
  (`github.com/drchristhorpe/histo_hmm`, branch `master`).
- Claude Code / Claude Desktop / Claude Science skill
  (`skills/histo-match/`) wrapping the CLI.
- Test suite (pytest) against real bundled reference sequences (no
  synthetic fixtures needed).
- `README.md`, `CLAUDE.md`, and design plan (`docs/PLAN.md`).
