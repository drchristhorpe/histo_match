---
name: histo-match
description: Match an MHC Class I protein sequence to known allele numbers (e.g. HLA-A*02:01), using profile HMM locus prediction to narrow the search space. Use when the user provides an amino-acid sequence and asks what allele/HLA type it is, wants it identified against known allele databases, or asks to resolve a sequence to an allele name/number.
---

# histo-match

`histo-match` is a CLI tool (from the `histo_match` package) that matches
an MHC Class I protein sequence against known allele numbers. Invoke it
with the Bash tool.

## When to use this skill

The user provides (or references) one or more amino-acid sequences and
asks:

- "what allele is this sequence?"
- "which HLA type does this correspond to?"
- to identify/resolve a sequence against known allele databases (IPD-IMGT/HLA, IPD-MHC, etc.)
- to batch-process a FASTA file of candidate sequences for allele assignment

This complements the [`histo-hmm`](https://github.com/drchristhorpe/histo_hmm)
skill, which only predicts the locus (e.g. `hla_a`) — use `histo-match`
when the user wants the specific **allele**, not just the locus.

## Checking availability

```bash
histo-match --help
```

If this fails with "command not found", install it first:

```bash
uv tool install .   # from a checkout of the histo_match repo
# or
pip install histo-match
```

(If working from a checkout of the `histo_match` source repo instead of an
installed package, use `uv run histo-match ...` there instead.)

## Usage

```bash
histo-match match --sequence "MAVMAPRTLVLLLSG..." [--plain] [--json-file PATH]
histo-match match --fasta sequences.fasta [--top 10] [--json-file PATH]
```

- `--sequence` — repeat for multiple inputs.
- `--fasta` — match every sequence in a FASTA file.
- `--top` — number of ranked candidate alleles to report when there's no
  exact match (default 5).
- `--candidate-loci` — number of HMM-ranked loci to search (default 3);
  raise this if a sequence is unexpectedly reported as no match — closely
  related loci (e.g. `BoLA-1`/`BoLA-2`) don't always rank the true locus in
  the default top 3.
- `--json-file PATH` — write full structured results to a JSON file.
- `--plain` — compact machine-friendly stdout instead of Rich panels/tables.

## Interpreting output

For each sequence, `histo-match` reports:

- **Exact match** vs **Approximate match** vs **No match**
- The predicted locus and its HMM confidence
- The best-assigned allele: `allele_slug` (e.g. `hla-a_02_01`), full
  `gene_allele_name`/`protein_allele_name` (e.g. `HLA-A*02:01`), and
  database `allele_id`, plus its identity (1.0 for an exact match)
- A ranked table of alternate candidate alleles

When a sequence is shared by multiple allele names (or several alleles tie
on approximate identity), `histo-match` always assigns the one with the
**lowest allele number** — report that assignment, but mention the tied
alternates from the "Matches" table if the user asks what else it could be.

## Example

```bash
$ histo-match match --sequence "GSHSLRYLHILVSRPGHGSDLYSSVGFLDDTQFVRFSSDA..." --plain
sequence_1: exact=True locus=sla_6 locus_confidence=1.0000 best=sla-6_02_01 identity=1.0000
```

Report the allele assignment back to the user in whatever form they asked
for (allele name, slug, a table across a FASTA batch, etc.) — this skill
only tells you how to obtain it.
