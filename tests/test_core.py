import json

import pytest

from histo_match import AlleleMatcher
from histo_match.reference_data import clean_sequence, default_data_dir

DATA_DIR = default_data_dir()

# -- sla_6.json : real locus with a single sequence shared by two allele
#    names (SLA-6*02:01 and SLA-6*06:01), and one the HMM classifies as
#    sla_6 at rank 1 with full confidence -- the natural test case for the
#    "assign the lowest allele number" rule on an exact match --
SLA_6_MULTI_ALLELE_SEQUENCE = (
    "GSHSLRYLHILVSRPGHGSDLYSSVGFLDDTQFVRFSSDAANPRVEPRAPWMEQEGREYWDRQTDIAKEHSKASRSNLR"
    "VIIGNHNHSQSESHSFLWVSGCDVGSDGRILRGYEQFSYDGDDYIVLNEDLRSWTAISTVAQIIRRKWEAEGVAEQYRAY"
    "LEIECVEWLRKYLEKGKDVLQRAVPPKTHVTRHPFYDNKVTLRCWALGFYPKEISLTWQRDGEDQTQDMELVETRPSGD"
    "GTFQKWAALVVPSGEEQSYTCQVQHEGLQEPLTLRWE"
)


@pytest.fixture(scope="module")
def matcher() -> AlleleMatcher:
    return AlleleMatcher()


def test_exact_match_assigns_the_lowest_allele_number(matcher: AlleleMatcher):
    result = matcher.match(SLA_6_MULTI_ALLELE_SEQUENCE)

    assert result.exact_match is True
    assert result.best_match is not None
    assert result.best_match.identity == 1.0
    assert result.best_match.gene_allele_name == "SLA-6*02:01"
    assert result.best_match.allele_slug == "sla-6_02_01"
    assert [m.gene_allele_name for m in result.matches] == [
        "SLA-6*02:01",
        "SLA-6*06:01",
    ]


def test_approximate_match_on_single_residue_mismatch(matcher: AlleleMatcher):
    mutated = list(SLA_6_MULTI_ALLELE_SEQUENCE)
    mutated[20] = "A" if mutated[20] != "A" else "G"
    mutated_sequence = "".join(mutated)

    result = matcher.match(mutated_sequence)

    assert result.exact_match is False
    assert result.best_match is not None
    assert 0.0 < result.best_match.identity < 1.0
    # Ranked matches must be sorted by descending identity.
    identities = [m.identity for m in result.matches]
    assert identities == sorted(identities, reverse=True)


def test_match_strips_whitespace_and_gap_characters(matcher: AlleleMatcher):
    padded = f"  {SLA_6_MULTI_ALLELE_SEQUENCE[:10]}-?{SLA_6_MULTI_ALLELE_SEQUENCE[10:]}  \n"

    result = matcher.match(padded)

    assert result.cleaned_sequence == clean_sequence(padded)
    assert result.exact_match is True


def test_match_batch_matches_individual_results(matcher: AlleleMatcher):
    hla_a = json.loads((DATA_DIR / "hla_a.json").read_text())
    sequences = [SLA_6_MULTI_ALLELE_SEQUENCE, clean_sequence(next(iter(hla_a)))]

    batch_results = matcher.match_batch(sequences)
    individual_results = [matcher.match(seq) for seq in sequences]

    assert [r.best_match for r in batch_results] == [r.best_match for r in individual_results]


def test_match_returns_no_match_without_crashing_on_nonsense_input(matcher: AlleleMatcher):
    result = matcher.match("AAAAAAAAAA")

    assert result.query_sequence == "AAAAAAAAAA"
    # Whatever the HMM makes of ten alanines, this must not raise.
    assert isinstance(result.exact_match, bool)
