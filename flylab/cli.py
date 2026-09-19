from __future__ import annotations

import json

import typer

from flylab.assays.taste import run_taste_assay
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


@app.command()
def occupancy(
    compound: str = typer.Argument(..., help="Library key, e.g. imidacloprid"),
    conc: float = typer.Option(..., "--conc", help="Concentration in mol/L"),
    json_out: bool = typer.Option(False, "--json"),
) -> None:
    result = compare_compound(compound, conc)
    if json_out:
        typer.echo(json.dumps(result, indent=2))
    else:
        _print_table(result)


@app.command()
def compare(
    compounds: list[str] = typer.Argument(...),
    conc: float = typer.Option(..., "--conc"),
) -> None:
    for name in compounds:
        _print_table(compare_compound(name, conc))
        typer.echo("-" * 64)


@app.command("list-drugs")
def list_drugs() -> None:
    lib = load_library()
    for key, spec in lib["compounds"].items():
        typer.echo(f"{key:16} {spec['name']}  ({spec.get('class', '')})")


@app.command()
def assay(
    compound: str = typer.Option("imidacloprid", "--compound"),
    conc: float = typer.Option(1e-6, "--conc"),
    sugar: float = typer.Option(150.0, "--sugar"),
    bitter: float = typer.Option(0.0, "--bitter"),
) -> None:
    nb = run_taste_assay(compound=compound, conc_M=conc, sugar_hz=sugar, bitter_hz=bitter)
    typer.echo(json.dumps(nb, indent=2))


@app.command()
def serve(host: str = "127.0.0.1", port: int = 8765) -> None:
    import uvicorn

    typer.echo(f"FlyLab bench → http://{host}:{port}")
    uvicorn.run("flylab.server:app", host=host, port=port, reload=False)


if __name__ == "__main__":
    app()
