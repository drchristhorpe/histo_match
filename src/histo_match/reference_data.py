from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

_GAP_CHARS = str.maketrans("", "", "-?")


def clean_sequence(sequence: str) -> str:
    """Uppercase and strip whitespace and `-`/`?` placeholder characters.

    Reference sequences in the bundled data use `-` to mark an
    unresolved/missing residue and `?` for an unknown one; `histo_hmm`
    strips the same characters before training/scoring, so query and
    reference sequences must be cleaned identically to compare like-for-like.
    """
    return "".join(sequence.split()).upper().translate(_GAP_CHARS)


def slugify_allele_name(name: str) -> str:
    """Derive a URL/filename-safe slug from an allele name.

    Lowercases the name and replaces `*`/`:` with `_`, e.g.
    "HLA-A*02:01" -> "hla-a_02_01", "Rano-E*u" -> "rano-e_u".
    """
    return name.lower().replace("*", "_").replace(":", "_")


def allele_sort_key(name: str) -> tuple:
    """Natural sort key for allele names so numeric fields compare numerically.

    Splits on runs of digits, e.g. "HLA-A*01:01:05" and "HLA-A*02:01" become
    comparable tuples; non-numeric nomenclature (e.g. "Rano-E*u") falls back
    to plain string comparison of its parts.
    """
    parts = re.split(r"(\d+)", name)
    return tuple(int(p) if p.isdigit() else p for p in parts)


def default_data_dir() -> Path:
    return Path(__file__).parent / "data" / "cytoplasmic_sequences"


class LocusReferenceData:
    """The cleaned exact-match lookup and ordered sequence list for one locus."""

    def __init__(self, locus: str, raw: dict) -> None:
        self.locus = locus

        # The source data's own "canonical_allele" field is sometimes an
        # empty dict (77 of 14,180 bundled entries) rather than missing
        # entirely, so it can't be trusted as-is. The canonical allele is
        # always self-derived here as the lowest-numbered entry in
        # "alleles", which is the one rule this must never get wrong.
        by_sequence: dict[str, dict] = {}
        for raw_sequence, entry in raw.items():
            cleaned = clean_sequence(raw_sequence)
            existing = by_sequence.get(cleaned)
            if existing is None:
                by_sequence[cleaned] = {"alleles": list(entry["alleles"])}
            else:
                # Two raw keys cleaned to the same sequence: merge allele lists.
                existing["alleles"].extend(
                    a for a in entry["alleles"] if a not in existing["alleles"]
                )

        for entry in by_sequence.values():
            entry["alleles"].sort(key=lambda a: allele_sort_key(a["gene_allele_name"]))
            entry["canonical_allele"] = entry["alleles"][0]

        self._by_sequence = by_sequence
        self.ordered_sequences = sorted(
            by_sequence,
            key=lambda seq: allele_sort_key(
                by_sequence[seq]["canonical_allele"]["gene_allele_name"]
            ),
        )

    def exact_lookup(self, cleaned_sequence: str) -> dict | None:
        return self._by_sequence.get(cleaned_sequence)

    def ordered_items(self):
        """Yield (sequence, entry) pairs ordered by ascending allele number."""
        for seq in self.ordered_sequences:
            yield seq, self._by_sequence[seq]


@lru_cache(maxsize=None)
def _load_locus_cached(data_dir: str, locus: str) -> LocusReferenceData:
    path = Path(data_dir) / f"{locus}.json"
    with path.open() as f:
        raw = json.load(f)
    return LocusReferenceData(locus, raw)


def load_locus(locus: str, data_dir: Path | str | None = None) -> LocusReferenceData:
    """Load (and cache) the reference data for one locus, e.g. "hla_a"."""
    resolved = Path(data_dir) if data_dir is not None else default_data_dir()
    return _load_locus_cached(str(resolved), locus)
