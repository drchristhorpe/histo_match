import json

from click.testing import CliRunner

from histo_match.cli import main

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


def test_match_requires_sequence_or_fasta():
    result = CliRunner().invoke(main, ["match"])

    assert result.exit_code != 0
    assert "Provide at least one --sequence or --fasta" in result.output


def test_match_sequence_plain_output():
    result = CliRunner().invoke(main, ["match", "--sequence", SLA_6_MULTI_ALLELE_SEQUENCE, "--plain"])

    assert result.exit_code == 0
    assert "exact=True" in result.output
    assert "sla-6_02_01" in result.output


def test_match_writes_json_file(tmp_path):
    json_path = tmp_path / "results.json"

    result = CliRunner().invoke(
        main,
        [
            "match",
            "--sequence",
            SLA_6_MULTI_ALLELE_SEQUENCE,
            "--plain",
            "--json-file",
            str(json_path),
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(json_path.read_text())
    assert len(payload) == 1
    assert payload[0]["exact_match"] is True
    assert payload[0]["best_match"]["allele_slug"] == "sla-6_02_01"


def test_match_fasta_input(tmp_path):
    fasta_path = tmp_path / "sequences.fasta"
    fasta_path.write_text(f">query_one\n{SLA_6_MULTI_ALLELE_SEQUENCE}\n")

    result = CliRunner().invoke(main, ["match", "--fasta", str(fasta_path), "--plain"])

    assert result.exit_code == 0
    assert "query_one: exact=True" in result.output
