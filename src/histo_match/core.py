from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import Levenshtein
from histo_hmm import MHCClassIClassifier

from .reference_data import allele_sort_key, clean_sequence, load_locus, slugify_allele_name


@dataclass(frozen=True)
class AlleleMatch:
    """One candidate allele assignment for a query sequence."""

    allele_slug: str
    gene_allele_name: str
    protein_allele_name: str
    allele_id: str
    locus: str
    source: str
    identity: float
    matched_sequence: str

    @classmethod
    def from_record(cls, record: dict, identity: float, matched_sequence: str) -> AlleleMatch:
        return cls(
            allele_slug=slugify_allele_name(record["gene_allele_name"]),
            gene_allele_name=record["gene_allele_name"],
            protein_allele_name=record["protein_allele_name"],
            allele_id=record["id"],
            locus=record["locus"],
            source=record["source"],
            identity=identity,
            matched_sequence=matched_sequence,
        )


@dataclass(frozen=True)
class MatchResult:
    """The outcome of matching one query sequence against the reference data."""

    query_sequence: str
    cleaned_sequence: str
    is_class_i: bool
    predicted_locus: str | None
    locus_confidence: float
    candidate_loci: list[tuple[str, float]]
    region_start: int
    region_end: int
    exact_match: bool
    best_match: AlleleMatch | None
    matches: list[AlleleMatch]


class AlleleMatcher:
    """Matches MHC Class I protein sequences to known allele numbers.

    Uses `histo_hmm` to narrow the search space to a handful of candidate
    loci, then matches against the reference sequences bundled for those
    loci only (exact match first, falling back to Levenshtein-ranked
    approximate matching).
    """

    def __init__(
        self,
        data_dir: str | Path | None = None,
        n_candidate_loci: int = 3,
        min_locus_probability: float = 0.0,
    ) -> None:
        self._classifier = MHCClassIClassifier()
        self._data_dir = Path(data_dir) if data_dir is not None else None
        self.n_candidate_loci = n_candidate_loci
        self.min_locus_probability = min_locus_probability

    def match(self, sequence: str, n_top: int = 5) -> MatchResult:
        cleaned = clean_sequence(sequence)
        classification = self._classifier.classify(
            cleaned, n_top=max(self.n_candidate_loci, 1), scan_constructs=True
        )
        trimmed = cleaned[classification.region_start : classification.region_end]

        candidate_loci = [
            (locus, probability)
            for locus, probability in classification.top_loci
            if probability >= self.min_locus_probability
        ][: self.n_candidate_loci]

        exact_match = False
        best_match: AlleleMatch | None = None
        matches: list[AlleleMatch] = []

        for locus, _probability in candidate_loci:
            reference = load_locus(locus, self._data_dir)
            entry = reference.exact_lookup(trimmed)
            if entry is not None:
                exact_match = True
                sorted_alleles = sorted(
                    entry["alleles"], key=lambda a: allele_sort_key(a["gene_allele_name"])
                )
                matches = [
                    AlleleMatch.from_record(allele, identity=1.0, matched_sequence=trimmed)
                    for allele in sorted_alleles
                ]
                best_match = matches[0]
                break

        if not exact_match:
            scored: list[tuple[float, dict, str]] = []
            for locus, _probability in candidate_loci:
                reference = load_locus(locus, self._data_dir)
                for ref_sequence, entry in reference.ordered_items():
                    identity = Levenshtein.ratio(trimmed, ref_sequence)
                    scored.append((identity, entry, ref_sequence))
            scored.sort(key=lambda item: -item[0])
            matches = [
                AlleleMatch.from_record(
                    item[1]["canonical_allele"], identity=item[0], matched_sequence=item[2]
                )
                for item in scored[:n_top]
            ]
            best_match = matches[0] if matches else None

        return MatchResult(
            query_sequence=sequence,
            cleaned_sequence=cleaned,
            is_class_i=classification.is_class_i,
            predicted_locus=classification.top_loci[0][0] if classification.top_loci else None,
            locus_confidence=classification.confidence,
            candidate_loci=candidate_loci,
            region_start=classification.region_start,
            region_end=classification.region_end,
            exact_match=exact_match,
            best_match=best_match,
            matches=matches,
        )

    def match_batch(self, sequences: list[str], n_top: int = 5) -> list[MatchResult]:
        return [self.match(sequence, n_top=n_top) for sequence in sequences]
