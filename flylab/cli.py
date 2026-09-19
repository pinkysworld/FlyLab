from __future__ import annotations
import json
from pathlib import Path
import typer
from flylab.assays.taste import run_taste_assay
from flylab.assays.wholens import run_wholens_assay
from flylab.assays.subgraph import run_subgraph_assay
from flylab.maps.malecns import download as download_malecns_files
from flylab.pharm.occupancy import compare_compound, load_library

app = typer.Typer(help="FlyLab virtual fly pharmacology bench")

def _print_table(result: dict) -> None:
    typer.echo(f"{result['compound']}  @  {result['concentration_M']:.2e} M")
    typer.echo(f"class: {result['class']}")
    typer.echo(f"{'receptor':28} {'occupancy':>10} {'EC50 (M)':>12}  direction")
    for row in result["receptors"]:
        typer.echo(f"{row['receptor']:28} {row['occupancy']:10.3f} {row['ec50_M']:12.2e}  {row['direction']}")
    typer.echo(result["disclaimer"])

@app.command()
def occupancy(compound: str, conc: float = typer.Option(..., "--conc"), json_out: bool = False):
    result = compare_compound(compound, conc)
    if json_out:
        typer.echo(json.dumps(result, indent=2))
    else:
        _print_table(result)

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

@app.command("extract-subgraph")
def extract_subgraph(
    hops: int = 1,
    min_weight: int = 5,
    types: str = typer.Option("MN9,DNp01", help="Comma-separated MaleCNS type names"),
):
    from flylab.maps.extract import extract_neighborhood
    type_list = tuple(t.strip() for t in types.split(",") if t.strip())
    payload = extract_neighborhood(
        types=type_list,
        hops=hops,
        min_weight=min_weight,
        out=Path("data/derived/malecns_named_neighborhood.json"),
    )
    typer.echo(f"seeds={payload['seeds']}")
    typer.echo(f"nodes={payload['n_nodes']} edges={payload['n_edges']} -> {payload['path']}")

@app.command("assay-cns")
def assay_cns(compound: str = "imidacloprid", conc: float = 1e-6):
    typer.echo(json.dumps(run_wholens_assay(compound, conc), indent=2, default=str))

@app.command("assay-subgraph")
def assay_subgraph(compound: str = "imidacloprid", conc: float = 1e-6):
    typer.echo(json.dumps(run_subgraph_assay(compound, conc), indent=2))

@app.command()
def serve(host: str = "127.0.0.1", port: int = 8765):
    import uvicorn
    typer.echo(f"FlyLab bench -> http://{host}:{port}")
    uvicorn.run("flylab.server:app", host=host, port=port, reload=False)

if __name__ == "__main__":
    app()
