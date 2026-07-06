import json
from pathlib import Path

from histo_match.reference_data import (
    LocusReferenceData,
    allele_sort_key,
    clean_sequence,
    default_data_dir,
    load_locus,
    slugify_allele_name,
)

DATA_DIR = default_data_dir()

# -- bola_1.json : real locus with a single sequence shared by two allele
#    names (BoLA-1*007:01:01:01 and :02), the natural test case for the
#    "assign the lowest allele number" rule --
BOLA_1_MULTI_ALLELE_SEQUENCE = (
    "GSHSLKYFHTAVSRPGDGEPRFITVGYVDDTQFVRFDSDAPDPRKEPRAPWIEKEGPEYWDRETRISKENTLVYRGSL"
    "NNLRGYYNQSEAGSHTFQQMYGCDVGPDGRLLRGFKQFAYDSRDYIALNEELRSWTAADTAAQITKRKWEAAGAAETWR"
    "NYLEGECVEWLRRYLENGKDTLLRADPPKAHVTHHPISDREVTLRCWALGFYPKEISLTWQRNGEDQTQDMELVETRPS"
    "GDGNFQKWAALVVPSGEEQRYTCHVQHEGLQEPLTLRWE"
)


def test_clean_sequence_strips_whitespace_case_and_gap_chars():
    assert clean_sequence(" mav-map?rt \n") == "MAVMAPRT"


def test_slugify_allele_name():
    assert slugify_allele_name("HLA-A*02:01") == "hla-a_02_01"
    assert slugify_allele_name("Rano-E*u") == "rano-e_u"
    assert slugify_allele_name("HLA-A*01:01:05") == "hla-a_01_01_05"


def test_allele_sort_key_compares_numeric_fields_numerically():
    # Plain string comparison would put "10" before "9"; the natural sort
    # key must not.
    assert allele_sort_key("HLA-A*2:9") < allele_sort_key("HLA-A*2:10")
    assert sorted(
        ["HLA-A*01:10", "HLA-A*01:02", "HLA-A*01:09"], key=allele_sort_key
    ) == ["HLA-A*01:02", "HLA-A*01:09", "HLA-A*01:10"]


def test_locus_reference_data_exact_lookup_picks_lowest_allele_number():
    raw = json.loads((DATA_DIR / "bola_1.json").read_text())
    reference = LocusReferenceData("bola_1", raw)

    entry = reference.exact_lookup(BOLA_1_MULTI_ALLELE_SEQUENCE)

    assert entry is not None
    names = {allele["gene_allele_name"] for allele in entry["alleles"]}
    assert names == {"BoLA-1*007:01:01:01", "BoLA-1*007:01:01:02"}
    assert entry["canonical_allele"]["gene_allele_name"] == "BoLA-1*007:01:01:01"


def test_locus_reference_data_merges_sequences_that_clean_to_the_same_string():
    raw = {
        "-MAVMAPRT": {
            "alleles": [{"gene_allele_name": "Test-A*02:01", "id": "X2", "locus": "A", "protein_allele_name": "Test-A*02:01", "source": "test"}],
            "canonical_allele": {"gene_allele_name": "Test-A*02:01", "id": "X2", "locus": "A", "protein_allele_name": "Test-A*02:01", "source": "test"},
        },
        "MAVMAPRT": {
            "alleles": [{"gene_allele_name": "Test-A*01:01", "id": "X1", "locus": "A", "protein_allele_name": "Test-A*01:01", "source": "test"}],
            "canonical_allele": {"gene_allele_name": "Test-A*01:01", "id": "X1", "locus": "A", "protein_allele_name": "Test-A*01:01", "source": "test"},
        },
    }

    reference = LocusReferenceData("test_locus", raw)

    entry = reference.exact_lookup("MAVMAPRT")
    assert entry is not None
    assert {a["gene_allele_name"] for a in entry["alleles"]} == {"Test-A*02:01", "Test-A*01:01"}
    assert entry["canonical_allele"]["gene_allele_name"] == "Test-A*01:01"


def test_ordered_items_ascend_by_allele_number():
    raw = {
        "SEQTWO": {
            "alleles": [{"gene_allele_name": "Test-A*02:01", "id": "X2", "locus": "A", "protein_allele_name": "Test-A*02:01", "source": "test"}],
            "canonical_allele": {"gene_allele_name": "Test-A*02:01", "id": "X2", "locus": "A", "protein_allele_name": "Test-A*02:01", "source": "test"},
        },
        "SEQONE": {
            "alleles": [{"gene_allele_name": "Test-A*01:01", "id": "X1", "locus": "A", "protein_allele_name": "Test-A*01:01", "source": "test"}],
            "canonical_allele": {"gene_allele_name": "Test-A*01:01", "id": "X1", "locus": "A", "protein_allele_name": "Test-A*01:01", "source": "test"},
        },
    }

    reference = LocusReferenceData("test_locus", raw)

    sequences = [seq for seq, _entry in reference.ordered_items()]
    assert sequences == ["SEQONE", "SEQTWO"]


def test_load_locus_caches_by_data_dir_and_locus():
    first = load_locus("bola_1")
    second = load_locus("bola_1")
    assert first is second
    assert isinstance(first, LocusReferenceData)
