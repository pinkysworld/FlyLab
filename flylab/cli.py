from __future__ import annotations

import json

import typer

from flylab.pharm.occupancy import compare_compound, load_library

app = typer.Typer(help="FlyLab — virtual fly pharmacology bench")


def _print_table(result: dict) -> None:
    typer.echo(f"{result['compound']}  @  {result['concentration_M']:.2e} M")
    typer.echo(f"class: {result['class']}")
    typer.echo("")
    typer.echo(f"{'receptor':28} {'occupancy':>10} {'EC50 (M)':>12}  direction")
    for row in result["receptors"]:
        typer.echo(
            f"{row['receptor']:28} {row['occupancy']:10.3f} {row['ec50_M']:12.2e}  {row['direction']}"
        )
    typer.echo("")
    typer.echo(result["disclaimer"])
    typer.echo("Sources:")
    for row in result["receptors"]:
        typer.echo(f"  - {row['receptor']}: {row['source']}")


@app.command()
def occupancy(
    compound: str = typer.Argument(..., help="Library key, e.g. imidacloprid"),
    conc: float = typer.Option(..., "--conc", help="Concentration in mol/L"),
    json_out: bool = typer.Option(False, "--json"),
) -> None:
    """Receptor occupancy. No connectome required."""
    result = compare_compound(compound, conc)
    if json_out:
        typer.echo(json.dumps(result, indent=2))
    else:
        _print_table(result)


@app.command()
def compare(
    compounds: list[str] = typer.Argument(..., help="Two or more library keys"),
    conc: float = typer.Option(..., "--conc", help="Concentration in mol/L"),
) -> None:
    """Side-by-side occupancy at one concentration."""
    for name in compounds:
        _print_table(compare_compound(name, conc))
        typer.echo("-" * 64)


@app.command("list-drugs")
def list_drugs() -> None:
    lib = load_library()
    for key, spec in lib["compounds"].items():
        typer.echo(f"{key:16} {spec['name']}  ({spec.get('class', '')})")


if __name__ == "__main__":
    app()
