"""Command line interface for histo_match."""

from __future__ import annotations

import json
from pathlib import Path

import click
from Bio import SeqIO
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from histo_match.core import AlleleMatch, AlleleMatcher, MatchResult

console = Console()


def _clean_for_output(name: str) -> str:
    return name.replace("\n", " ").strip() or "sequence"


def _format_matches(matches: list[AlleleMatch]) -> Table:
    table = Table(title="Matches", show_header=True, header_style="bold cyan")
    table.add_column("Rank", justify="right")
    table.add_column("Allele slug", style="green")
    table.add_column("Allele name")
    table.add_column("Locus")
    table.add_column("Identity", justify="right")
    for index, match in enumerate(matches, start=1):
        table.add_row(
            str(index), match.allele_slug, match.gene_allele_name, match.locus, f"{match.identity:.4f}"
        )
    return table


def _render_result(name: str, result: MatchResult) -> None:
    if result.best_match is not None:
        status = "Exact match" if result.exact_match else "Approximate match"
        color = "green" if result.exact_match else "yellow"
    else:
        status = "No match"
        color = "red"

    summary = Table.grid(padding=(0, 2))
    summary.add_column(style="bold")
    summary.add_column()
    summary.add_row("Sequence", name)
    summary.add_row("Call", f"[{color}]{status}[/{color}]")
    summary.add_row("Predicted locus", result.predicted_locus or "none")
    summary.add_row("Locus confidence", f"{result.locus_confidence:.4f}")
    summary.add_row("Region", f"{result.region_start}-{result.region_end}")
    if result.best_match is not None:
        summary.add_row("Best allele", result.best_match.allele_slug)
        summary.add_row("Identity", f"{result.best_match.identity:.4f}")

    console.print(Panel(summary, title=name, expand=False))
    if result.matches:
        console.print(_format_matches(result.matches))


def _render_plain_result(name: str, result: MatchResult) -> None:
    best = result.best_match
    best_desc = f"{best.allele_slug} identity={best.identity:.4f}" if best is not None else "none"
    click.echo(
        f"{name}: exact={bool(result.exact_match)} locus={result.predicted_locus} "
        f"locus_confidence={float(result.locus_confidence):.4f} best={best_desc}"
    )


def _allele_match_to_dict(match: AlleleMatch) -> dict[str, object]:
    return {
        "allele_slug": match.allele_slug,
        "gene_allele_name": match.gene_allele_name,
        "protein_allele_name": match.protein_allele_name,
        "allele_id": match.allele_id,
        "locus": match.locus,
        "source": match.source,
        "identity": float(match.identity),
    }


def _result_to_dict(name: str, result: MatchResult) -> dict[str, object]:
    return {
        "name": name,
        "is_class_i": bool(result.is_class_i),
        "predicted_locus": result.predicted_locus,
        "locus_confidence": float(result.locus_confidence),
        "candidate_loci": [
            {"locus": locus, "probability": float(probability)}
            for locus, probability in result.candidate_loci
        ],
        "region_start": int(result.region_start),
        "region_end": int(result.region_end),
        "exact_match": bool(result.exact_match),
        "best_match": _allele_match_to_dict(result.best_match) if result.best_match else None,
        "matches": [_allele_match_to_dict(match) for match in result.matches],
    }


@click.group()
def main() -> None:
    """Match MHC Class I allele sequences to known allele numbers."""


@main.command("match")
@click.option(
    "--sequence",
    "sequences",
    multiple=True,
    help="Amino-acid sequence to match. Repeat for multiple inputs.",
)
@click.option(
    "--fasta",
    "fasta_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Match every sequence in a FASTA file.",
)
@click.option(
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory containing per-locus reference JSON (defaults to the bundled data).",
)
@click.option(
    "--candidate-loci",
    "n_candidate_loci",
    type=int,
    default=3,
    show_default=True,
    help="Number of HMM-ranked loci to search.",
)
@click.option("--min-locus-probability", type=float, default=0.0, show_default=True)
@click.option(
    "--top",
    "n_top",
    type=int,
    default=5,
    show_default=True,
    help="Number of ranked allele matches to report when there is no exact match.",
)
@click.option(
    "--plain",
    is_flag=True,
    help="Print plain text output instead of Rich-formatted panels and tables.",
)
@click.option(
    "--json-file",
    "json_output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional JSON filename for full match results.",
)
def match_cmd(
    sequences: tuple[str, ...],
    fasta_path: Path | None,
    data_dir: Path | None,
    n_candidate_loci: int,
    min_locus_probability: float,
    n_top: int,
    plain: bool,
    json_output: Path | None,
) -> None:
    """Match one or more sequences to known allele numbers."""
    if not sequences and fasta_path is None:
        raise click.UsageError("Provide at least one --sequence or --fasta")

    inputs: list[tuple[str, str]] = []
    for i, seq in enumerate(sequences, start=1):
        inputs.append((f"sequence_{i}", seq))

    if fasta_path is not None:
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            inputs.append((_clean_for_output(record.id), str(record.seq)))

    matcher = AlleleMatcher(
        data_dir=data_dir,
        n_candidate_loci=n_candidate_loci,
        min_locus_probability=min_locus_probability,
    )

    payload: list[dict[str, object]] = []
    for name, seq in inputs:
        result = matcher.match(seq, n_top=n_top)
        if plain:
            _render_plain_result(name, result)
        else:
            _render_result(name, result)
        payload.append(_result_to_dict(name, result))

    if json_output is not None:
        json_output.parent.mkdir(parents=True, exist_ok=True)
        with open(json_output, "w") as f:
            json.dump(payload, f, indent=2)
        console.print(f"[bold green]Wrote JSON results to[/bold green] {json_output}")


if __name__ == "__main__":
    main()
