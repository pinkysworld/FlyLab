from __future__ import annotations
import json
import typer
from flylab.assays.taste import run_taste_assay
from flylab.assays.wholens import run_wholens_assay
from flylab.maps.malecns import download as download_malecns_files
from flylab.pharm.occupancy import compare_compound, load_library

app = typer.Typer(help="FlyLab — virtual fly pharmacology bench")

def _print_table(result: dict) -> None:
    typer.echo(f"{result['compound']}  @  {result['concentration_M']:.2e} M")
    typer.echo(f"class: {result['class']}\n")
    typer.echo(f"{'receptor':28} {'occupancy':>10} {'EC50 (M)':>12}  direction")
    for row in result["receptors"]:
        typer.echo(f"{row['receptor']:28} {row['occupancy']:10.3f} {row['ec50_M']:12.2e}  {row['direction']}")
    typer.echo("\n" + result["disclaimer"])

@app.command()
def occupancy(compound: str, conc: float = typer.Option(..., "--conc"), json_out: bool = typer.Option(False, "--json")):
    result = compare_compound(compound, conc)
    typer.echo(json.dumps(result, indent=2) if json_out else None) if json_out else _print_table(result)
    if json_out:
        typer.echo(json.dumps(compare_compound(compound, conc), indent=2))

@app.command()
def compare(compounds: list[str], conc: float = typer.Option(..., "--conc")):
    for name in compounds:
        _print_table(compare_compound(name, conc))

@app.command("list-drugs")
def list_drugs():
    for key, spec in load_library()["compounds"].items():
        typer.echo(f"{key:16} {spec['name']}")

@app.command()
def assay(compound: str = "imidacloprid", conc: float = 1e-6, sugar: float = 150.0, bitter: float = 0.0):
    typer.echo(json.dumps(run_taste_assay(compound, conc, sugar, bitter), indent=2))

@app.command("download-malecns")
def download_malecns(full: bool = typer.Option(False, "--full")):
    typer.echo(str(download_malecns_files(kind="full" if full else "atlas")))

@app.command("assay-cns")
def assay_cns(compound: str = "imidacloprid", conc: float = 1e-6):
    typer.echo(json.dumps(run_wholens_assay(compound, conc), indent=2, default=str))

@app.command()
def serve(host: str = "127.0.0.1", port: int = 8765):
    import uvicorn
    typer.echo(f"FlyLab bench → http://{host}:{port}")
    uvicorn.run("flylab.server:app", host=host, port=port, reload=False)

if __name__ == "__main__":
    app()
