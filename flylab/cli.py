"""FlyLab command line.

Every command is a thin wrapper over the same library functions the HTTP API
uses, so a result obtained at the terminal and one obtained in the bench UI are
the same numbers.  Commands that need a module another agent is still writing
import it lazily and print a friendly message instead of a traceback.

``--json`` is available on every command that prints a table, and emits exactly
the payload the matching endpoint returns.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import typer
from click import exceptions as click_exceptions

from flylab.assays.subgraph import run_subgraph_assay
from flylab.assays.taste import run_taste_assay
from flylab.assays.wholens import run_wholens_assay
from flylab.maps.malecns import download as download_malecns_files
from flylab.pharm.occupancy import compare_compound, load_library

app = typer.Typer(help="FlyLab virtual fly pharmacology bench", no_args_is_help=True)
experiment_app = typer.Typer(help="Batch experiment designs", no_args_is_help=True)
graph_app = typer.Typer(help="Neighborhood graph inspection", no_args_is_help=True)
app.add_typer(experiment_app, name="experiment")
app.add_typer(graph_app, name="graph")

JSON_OPT = typer.Option(False, "--json", help="Print the raw JSON payload instead of a table.")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _dump(payload: Any) -> None:
    typer.echo(json.dumps(payload, indent=2, default=str))


def _fmt(value: Any, width: int = 10, digits: int = 3) -> str:
    if value is None:
        return "-".rjust(width)
    if isinstance(value, float):
        return f"{value:{width}.{digits}f}"
    return str(value).rjust(width)


def _sci(value: Any, width: int = 10, digits: int = 2) -> str:
    """Scientific notation, right-aligned; a missing value prints 'not modelled'."""
    if value is None:
        return "not modelled".rjust(width)
    try:
        return f"{float(value):{width}.{digits}e}"
    except (TypeError, ValueError):
        return str(value).rjust(width)


def _warnings(payload: dict[str, Any]) -> None:
    for w in payload.get("warnings") or []:
        typer.secho(f"  ! {w}", fg=typer.colors.YELLOW)


def _print_table(result: dict) -> None:
    typer.echo(f"{result['compound']}  @  {result['concentration_M']:.2e} M")
    typer.echo(f"class: {result['class']}")
    typer.echo(f"{'receptor':28} {'engagement':>10} {'value (M)':>12} {'type':>6}  direction")
    for row in result["receptors"]:
        # schema v3: a not-modelled row prints N/A, never 0.000
        eng = row.get("engagement", row.get("occupancy"))
        value = row.get("param_value_M", row.get("ec50_M"))
        eng_s = "       n/a" if eng is None else f"{eng:10.3f}"
        val_s = "not modelled" if value is None else f"{value:12.2e}"
        typer.echo(
            f"{row['receptor']:28} {eng_s} {val_s} {str(row.get('param_type', '')):>6}  {row['direction']}"
        )
    if result.get("engagement_is_not_occupancy"):
        typer.echo(result["engagement_is_not_occupancy"])
    typer.echo(result["disclaimer"])


def _readout_table(nb: dict[str, Any], keys: tuple[str, ...]) -> None:
    r = nb.get("readouts", {})
    typer.echo(f"{nb.get('assay')}  compound={nb.get('compound')}  conc={nb.get('concentration_M')}")
    for k in keys:
        if k in r:
            typer.echo(f"  {k:24} {_fmt(r[k])}")
    gains = nb.get("gains") or {}
    if gains:
        typer.echo("  gains " + "  ".join(f"{k}={v:.3f}" for k, v in gains.items()))
    _warnings(nb)


def _bridge(route: str, payload: dict[str, Any]) -> Any:
    """Run one bench route through the transport-free bridge.

    The CLI, the HTTP API and the browser build all end up in the same
    function, so a number obtained at the terminal is the number the dashboard
    shows.
    """
    from flylab.browser import bridge

    out = bridge.call(route, payload)
    if isinstance(out, dict) and isinstance(out.get("error"), dict):
        typer.secho(
            f"{route}: {out['error'].get('detail')}", fg=typer.colors.YELLOW, err=True
        )
        raise typer.Exit(code=2)
    return out


def _estimate_line(payload: dict[str, Any]) -> None:
    est = payload.get("runtime_estimate") or {}
    if est.get("estimate_s") is not None:
        typer.secho(
            f"  ~{est['estimate_s']} s estimated"
            + (f" - {est['note']}" if est.get("note") else ""),
            fg=typer.colors.BLUE,
        )


def _lazy_call(dotted: str, name: str, *args: Any, **kw: Any) -> Any:
    """Import and call, or exit(2) with a message naming the missing module."""
    try:
        module = __import__(dotted, fromlist=[name])
        fn = getattr(module, name)
    except (ImportError, AttributeError):
        typer.secho(
            f"{dotted}.{name} is not available in this build yet "
            "(it ships with the v0.5 analysis layer).",
            fg=typer.colors.YELLOW,
        )
        raise typer.Exit(code=2)
    return fn(*args, **kw)


# --------------------------------------------------------------------------
# library / occupancy
# --------------------------------------------------------------------------
@app.command()
def occupancy(
    compound: str,
    conc: float = typer.Option(..., "--conc", help="Concentration in molar."),
    json_out: bool = JSON_OPT,
):
    """Insect vs vertebrate receptor occupancy for one compound at one dose."""
    result = compare_compound(compound, conc)
    _dump(result) if json_out else _print_table(result)


@app.command()
def compare(
    compounds: list[str],
    conc: float = typer.Option(..., "--conc"),
    dependence: bool = typer.Option(
        True, "--dependence/--no-dependence", help="Include the topology-dependence column."
    ),
    n: int = typer.Option(20, "--n", help="Permutation shuffles for the dependence column."),
    occupancy_tables: bool = typer.Option(
        True, "--occupancy/--no-occupancy", help="Also print the per-compound occupancy tables."
    ),
    json_out: bool = JSON_OPT,
):
    """Compare compounds: receptor selectivity beside circuit selectivity.

    The same payload the dashboard's compare screen uses. It makes visible that
    the highest receptor selectivity is not the highest circuit selectivity.
    """
    out = _bridge(
        "/api/compare",
        {
            "compounds": list(compounds),
            "conc_M": conc,
            "include_dependence": dependence,
            "n": n,
        },
    )
    if json_out:
        _dump(out)
        return
    est = out.get("runtime_estimate") or {}
    typer.echo(f"compare @ {conc:.2e} M  (estimated {est.get('estimate_s')} s)")
    head = (
        f"{'compound':16} {'target':22} {'insect':>7} {'vert':>7} "
        f"{'recSI':>7} {'cirSI':>7} {'gap':>7} {'circuit d%':>11}  topology"
    )
    typer.echo(head)
    for row in out.get("rows") or []:
        typer.echo(
            f"{str(row.get('name'))[:16]:16} "
            f"{str(row.get('target_receptor') or '-')[:22]:22} "
            f"{_fmt(row.get('insect_engagement'), 7, 3)} "
            f"{_fmt(row.get('vertebrate_engagement'), 7, 3)} "
            f"{_fmt(row.get('receptor_si_log10'), 7, 2)} "
            f"{_fmt(row.get('circuit_si_log10'), 7, 2)} "
            f"{_fmt(row.get('si_gap_circuit_minus_receptor'), 7, 2)} "
            f"{_fmt(row.get('circuit_delta_percent'), 11, 1)}  "
            f"{row.get('topology_dependence') or '-'}"
        )
    for line in out.get("findings") or []:
        typer.secho("  * " + line, fg=typer.colors.CYAN)
    _warnings(out)
    if occupancy_tables:
        for name in compounds:
            typer.echo("")
            _print_table(compare_compound(name, conc))


@app.command("list-drugs")
def list_drugs(json_out: bool = JSON_OPT):
    """List every compound key in the library."""
    lib = load_library()["compounds"]
    if json_out:
        _dump([{"key": k, "name": v["name"], "class": v.get("class")} for k, v in lib.items()])
        return
    for key, spec in lib.items():
        typer.echo(f"{key:16} {spec['name']:28} {spec.get('class') or ''}")


@app.command()
def meta(json_out: bool = JSON_OPT):
    """Bench metadata: version, maps, library hash, graphs, mechanism rules."""
    from flylab.server import meta as meta_payload

    payload = meta_payload()
    if json_out:
        _dump(payload)
        return
    typer.echo(f"FlyLab {payload['version']}  notebook schema {payload['notebook_version']}")
    typer.echo(f"map   {payload['map']['id']}")
    lib = payload["library"]
    typer.echo(f"library {lib['version']}  sha256 {lib['sha256'][:16]}...  {lib['n_compounds']} compounds")
    for name, g in payload["graphs"].items():
        if not g.get("available"):
            typer.secho(f"graph {name:12} MISSING ({g.get('error')})", fg=typer.colors.YELLOW)
            continue
        seeds = ", ".join(f"{t}={n}" for t, n in g["seeds"].items())
        typer.echo(f"graph {name:12} {g['n_nodes']:6} nodes {g['n_edges']:7} edges  seeds: {seeds}")
    typer.echo(f"mechanism rules: {len(payload['mechanisms'])}")


# --------------------------------------------------------------------------
# assays
# --------------------------------------------------------------------------
@app.command()
def assay(
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    sugar: float = 150.0,
    bitter: float = 0.0,
    json_out: bool = JSON_OPT,
):
    """Reduced taste -> MN9 assay (reduced_taste_v0)."""
    nb = run_taste_assay(compound, conc, sugar, bitter)
    if json_out:
        _dump(nb)
    else:
        _readout_table(
            nb, ("g_ach", "mn9_sugar_hz", "mn9_sugar_bitter_hz", "mn9_vehicle_sugar_hz", "bitter_veto_ratio")
        )


@app.command("assay-cns")
def assay_cns(compound: str = "imidacloprid", conc: float = 1e-6, json_out: bool = JSON_OPT):
    """Whole-CNS census assay on the MaleCNS transmitter counts."""
    nb = run_wholens_assay(compound, conc)
    if json_out:
        _dump(nb)
        return
    r = nb["readouts"]
    typer.echo(f"traced cells: {r['n_traced']}  weights_present={r['weights_present']}")
    typer.echo(f"{'index':28} {'vehicle':>10} {'treated':>10}")
    for k in ("cns_excitation_index", "cns_inhibition_index", "excitation_inhibition_ratio"):
        typer.echo(f"{k:28} {_fmt(r['vehicle'][k])} {_fmt(r['treated'][k])}")
    _warnings(nb)


@app.command("assay-subgraph")
def assay_subgraph(
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    graph: str = typer.Option("named", "--graph", help="named | taste_motor | path to JSON"),
    drive_hz: float = 40.0,
    steps: int = 80,
    json_out: bool = JSON_OPT,
):
    """Rate-model assay on a committed MaleCNS neighborhood."""
    nb = run_subgraph_assay(compound, conc, drive_hz=drive_hz, steps=steps, graph=graph)
    if json_out:
        _dump(nb)
    else:
        _readout_table(nb, ("n_nodes", "n_edges", "mean_hz", "max_hz", "mn9_hz", "dnp01_hz"))


@app.command("assay-spiking")
def assay_spiking(
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    graph: str = typer.Option("named", "--graph"),
    drive_hz: float = 40.0,
    t_ms: float = 500.0,
    seed: int = 0,
    json_out: bool = JSON_OPT,
):
    """Shiu-style LIF spiking assay on a committed neighborhood."""
    from flylab.assays.spiking import run_spiking_assay

    nb = run_spiking_assay(compound, conc, drive_hz=drive_hz, t_ms=t_ms, seed=seed, graph=graph)
    if json_out:
        _dump(nb)
        return
    r = nb["readouts"]
    typer.echo(f"LIF {r['t_ms']:.0f} ms, {r['n_nodes']} cells, {r['n_spikes']} spikes")
    typer.echo(f"{'readout':24} {'vehicle':>10} {'treated':>10}")
    for k in ("mn9_hz", "dnp01_hz", "mean_hz"):
        typer.echo(f"{k:24} {_fmt(r['vehicle'].get(k))} {_fmt(r.get(k))}")
    _warnings(nb)


@app.command("assay-taste-map")
def assay_taste_map(
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    sugar: float = 150.0,
    bitter: float = 0.0,
    engine: str = typer.Option("rate", "--engine", help="rate | lif"),
    seed: int = 0,
    json_out: bool = JSON_OPT,
):
    """Map-extracted labellar GRN -> MN9 assay on the taste_motor graph."""
    from flylab.assays.taste_map import run_taste_map_assay

    nb = run_taste_map_assay(compound, conc, sugar_hz=sugar, bitter_hz=bitter, engine=engine, seed=seed)
    if json_out:
        _dump(nb)
    else:
        _readout_table(
            nb,
            (
                "engine",
                "n_sweet_grn",
                "n_bitter_grn",
                "mn9_sugar_hz",
                "mn9_sugar_bitter_hz",
                "mn9_vehicle_sugar_hz",
                "mn9_no_drive_hz",
                "bitter_veto_ratio",
            ),
        )


# --------------------------------------------------------------------------
# uncertainty / analysis
# --------------------------------------------------------------------------
@app.command()
def ensemble(
    assay: str = typer.Option("subgraph", "--assay", help="subgraph | spiking | taste | taste_map"),
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    n_rep: int = 8,
    ec50_sd_log10: float = 0.3,
    seed: int = 0,
    json_out: bool = JSON_OPT,
):
    """Replicate an assay and report a Monte-Carlo credible interval."""
    from flylab.assays.ensemble import run_ensemble

    nb = run_ensemble(assay, compound, conc, n_rep=n_rep, ec50_sd_log10=ec50_sd_log10, seed=seed)
    if json_out:
        _dump(nb)
        return
    u = nb["uncertainty"]
    typer.echo(f"{assay}  {compound} @ {conc:.2e} M  n_rep={u['n_rep']}")
    typer.echo(f"{'readout':16} {'mean':>10} {'sd':>10} {'lo':>10} {'hi':>10}")
    for k, (lo, hi) in u["ci"].items():
        typer.echo(f"{k:16} {_fmt(u['mean'].get(k))} {_fmt(u['sd'].get(k))} {_fmt(lo)} {_fmt(hi)}")
    _warnings(nb)


@app.command()
def sensitivity(
    assay: str = typer.Option("subgraph", "--assay"),
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    readout: str = "mean_hz",
    factor: float = 2.0,
    json_out: bool = JSON_OPT,
):
    """One-at-a-time tornado: which parameter moves the readout most."""
    from flylab.assays.ensemble import sensitivity as sensitivity_fn

    rows = sensitivity_fn(assay, compound, conc, factor=factor, readout=readout)
    if json_out:
        _dump(rows)
        return
    typer.echo(f"{assay}  {compound} @ {conc:.2e} M  readout={readout}  x/{factor}")
    typer.echo(f"{'parameter':20} {'low':>10} {'base':>10} {'high':>10} {'span':>10}")
    for row in rows:
        typer.echo(
            f"{str(row.get('param')):20} {_fmt(row.get('low'))} {_fmt(row.get('base'))} "
            f"{_fmt(row.get('high'))} {_fmt(row.get('span'))}"
        )


@app.command()
def ic50(
    assay: str = typer.Option("subgraph", "--assay"),
    compound: str = "imidacloprid",
    readout: str = "mean_hz",
    n_boot: int = 200,
    n_rep: int = 4,
    seed: int = 0,
    json_out: bool = JSON_OPT,
):
    """Hill fit of the model's own dose-response, with a bootstrap CI."""
    from flylab.assays.ensemble import circuit_ic50

    out = circuit_ic50(assay, compound, readout=readout, n_boot=n_boot, n_rep=n_rep, seed=seed)
    if json_out:
        _dump(out)
        return
    fit = out.get("fit") or {}
    typer.echo(f"{assay}  {compound}  readout={readout}")
    if fit.get("error"):
        typer.secho(f"  fit failed: {fit['error']}", fg=typer.colors.RED)
    else:
        typer.echo(f"  model IC50 {fit.get('ic50'):.3e} M   slope {_fmt(fit.get('slope'))}  r2 {_fmt(fit.get('r2'))}")
        ci = (out.get("ci") or {}).get("ic50")
        if ci:
            typer.echo(f"  bootstrap CI  {ci[0]:.3e} .. {ci[1]:.3e} M  ({out['bootstrap'].get('n_ok')} ok)")
    _warnings(out)


@app.command()
def exposure(
    compound: str = "imidacloprid",
    dose: float = typer.Option(1.0, "--dose", help="nanomoles delivered to one fly"),
    route: str = typer.Option("feeding", "--route", help="feeding | topical | bath"),
    t_h: float = 24.0,
    json_out: bool = JSON_OPT,
):
    """One-compartment exposure: C(t), AUC, Cmax, Tmax and occupancy(t)."""
    from flylab.pharm.exposure import exposure_profile

    out = exposure_profile(compound, dose, route, t_h=t_h)
    if json_out:
        _dump(out)
        return
    typer.echo(f"{out['compound']}  {out['route']}  {out['dose_nmol']} nmol")
    typer.echo(f"  Cmax {out['cmax']:.3e} M at Tmax {out['tmax']:.2f} h   AUC {out['auc']:.3e} M*h")
    _warnings(out)


@app.command("null-panel")
def null_panel(
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    n: int = typer.Option(8, "--n", help="Shuffles per mode. Each one re-runs the circuit."),
    seed: int = 0,
    modes: str = typer.Option(
        "", "--modes", help="Comma-separated subset of the shuffling modes (default: all)."
    ),
    json_out: bool = JSON_OPT,
):
    """Null-model panel: real effect against shuffled-network nulls."""
    kw = {"compound": compound, "conc_M": conc, "n": n, "seed": seed}
    if modes.strip():
        kw["modes"] = tuple(m.strip() for m in modes.split(",") if m.strip())
    out = _lazy_call("flylab.analysis.nullmodels", "null_panel", **kw)
    if json_out:
        _dump(out)
        return
    typer.echo(f"{compound} @ {conc:.2e} M  n={out.get('n')}")
    typer.echo(f"{'mode':28} {'real':>10} {'null mean':>10} {'null sd':>10} {'z':>8} {'p':>8}")
    for row in out.get("rows", []):
        typer.echo(
            f"{str(row.get('mode')):28} {_fmt(row.get('real_effect'))} {_fmt(row.get('null_mean'))} "
            f"{_fmt(row.get('null_sd'))} {_fmt(row.get('z'), 8, 2)} {_fmt(row.get('p_two_sided'), 8, 4)}"
        )
    _warnings(out)


@app.command()
def selectivity(conc: float = 1e-6, json_out: bool = JSON_OPT):
    """Receptor selectivity table across the library at one dose."""
    out = _lazy_call("flylab.analysis.selectivity", "receptor_selectivity_table", conc_M=conc)
    if json_out:
        _dump(out)
        return
    rows = out.get("rows", out) if isinstance(out, dict) else out
    typer.echo(f"{'compound':20} {'pair':10} {'insect occ':>11} {'vert occ':>10} {'SI log10':>9}")
    for row in rows if isinstance(rows, list) else []:
        typer.echo(
            f"{str(row.get('compound') or row.get('key')):20} {str(row.get('receptor_pair') or ''):10} "
            f"{_fmt(row.get('insect_occupancy'), 11)} {_fmt(row.get('vertebrate_occupancy'))} "
            f"{_fmt(row.get('receptor_si_log10'), 9, 2)}"
        )
    if isinstance(out, dict):
        _warnings(out)


@app.command()
def predictions(n_rep: int = 4, seed: int = 0, json_out: bool = JSON_OPT):
    """Falsifiable predictions table from the v0.5 analysis layer."""
    out = _lazy_call("flylab.analysis.predictions", "prediction_table", n_rep=n_rep, seed=seed)
    if json_out:
        _dump(out)
        return
    rows = out.get("rows", []) if isinstance(out, dict) else out
    for row in rows:
        typer.echo(f"{row.get('id', '?')}: {row.get('statement') or row.get('hypothesis') or row}")
    if isinstance(out, dict):
        _warnings(out)


# --------------------------------------------------------------------------
# dashboard / decision layer
# --------------------------------------------------------------------------
@app.command()
def dashboard(
    compound: str = typer.Argument("imidacloprid"),
    conc: float = typer.Option(1e-6, "--conc"),
    graph: str = typer.Option("named", "--graph"),
    dependence: bool = typer.Option(True, "--dependence/--no-dependence"),
    n: int = typer.Option(20, "--n", help="Permutation shuffles for the wiring verdict."),
    json_out: bool = JSON_OPT,
):
    """The decision dashboard for one compound at one dose, as one payload."""
    out = _bridge(
        "/api/dashboard",
        {
            "compound": compound,
            "conc_M": conc,
            "graph": graph,
            "include_dependence": dependence,
            "n": n,
        },
    )
    if json_out:
        _dump(out)
        return
    c = out["compound"]
    typer.secho(f"{c['name']}  -  {conc:.2e} M", bold=True)
    typer.echo(f"  class {c.get('class') or '-'}   target {c.get('target_receptor') or '-'}  ({c.get('mode') or '-'})")
    h = out["headline"]
    cov = out["coverage"]
    typer.echo(
        f"  insect engagement {_fmt(h['insect_engagement']['value'], 6, 3)} "
        f"[{h['insect_engagement']['classification']}] at {h['insect_engagement']['receptor']}"
        f"  ({h['insect_engagement'].get('engagement_model')}, "
        f"{h['insect_engagement'].get('evidence_distance_label')})"
    )
    typer.echo(
        f"  vertebrate engagement {_fmt(h['vertebrate_engagement']['value'], 6, 3)} "
        f"[{h['vertebrate_engagement']['classification']}] at {h['vertebrate_engagement']['receptor']}"
        f"  ({h['vertebrate_engagement'].get('engagement_model')}, "
        f"{h['vertebrate_engagement'].get('evidence_distance_label')})"
    )
    rs = h["receptor_selectivity"]
    typer.echo(
        f"  receptor selectivity {_fmt(rs.get('ratio_vert_over_insect'), 8, 1)}x "
        f"({rs.get('pair')})   evidence {h['evidence_tier']['value']} "
        f"({cov['n_sourced']} sourced / {cov['n_not_modelled']} not modelled)"
    )
    sel = out["selectivity"]
    typer.echo("\nselectivity (values, not a verdict)")
    typer.echo(
        f"  insect {_fmt(sel.get('insect_engagement'), 6, 3)} vs vertebrate "
        f"{_fmt(sel.get('vertebrate_engagement'), 6, 3)}  "
        f"difference {_fmt(sel.get('engagement_difference'), 6, 3)}"
    )
    lim = sel.get("limiting_vertebrate_receptor") or {}
    typer.echo(
        f"  vertebrate reaches {sel['occ_limit']:.0%} engagement at "
        f"{_sci(sel.get('vertebrate_limit_conc_M'))} M"
        + (f"; limiting receptor {lim.get('receptor')}" if lim else "")
    )
    typer.echo("\ncircuit")
    for key, r in (out["circuit"].get("readouts") or {}).items():
        typer.echo(
            f"  {r['label']:20} {r['vehicle']:8.3f} -> {r['treated']:8.3f} Hz "
            f"({r['percent']:+.1f}%) {r['direction']}"
        )
    dep = out.get("dependence")
    if dep:
        typer.echo(f"  {dep['verdict']}")
        for m in dep["modes"]:
            floor = " (at resolution floor)" if m["at_resolution_floor"] else ""
            typer.echo(
                f"    {m['mode']:26} p={m['p_two_sided']:.4f} "
                f"(resolution {m['p_resolution']:.4f}){floor} "
                f"{'beaten' if m['beats_null'] else 'not beaten'}"
            )
    typer.echo("\nwhat can I trust?")
    for row in out["trust"]["rows"]:
        typer.echo(f"  {row['area']:34} {row['status']:34} [{row['classification']}]")
    typer.echo(f"  {out['trust']['note']}")
    typer.echo("\nwhy this happened")
    typer.echo("  " + out["why"]["text"])
    _warnings(out)


@app.command()
def dependence(
    compound: str = typer.Argument("imidacloprid"),
    conc: float = typer.Option(1e-6, "--conc"),
    n: int = typer.Option(20, "--n"),
    graph: str | None = typer.Option(None, "--graph"),
    landscape: bool = typer.Option(False, "--landscape", help="The compound x conc landscape."),
    run: bool = typer.Option(False, "--run", help="Actually run the landscape (it is slow)."),
    json_out: bool = JSON_OPT,
):
    """Is this effect a wiring result? Permutation nulls over the MaleCNS cut."""
    if landscape:
        out = _bridge(
            "/api/dependence/landscape",
            {"n": n, "graph": graph, "estimate_only": not run},
        )
        if json_out:
            _dump(out)
            return
        _estimate_line(out)
        if out.get("estimate_only"):
            typer.echo("  pass --run to compute it")
            return
        _dump(out)
        return
    out = _bridge(
        "/api/dependence",
        {"compound": compound, "conc_M": conc, "n": n, "graph": graph},
    )
    if json_out:
        _dump(out)
        return
    _estimate_line(out)
    typer.secho(f"{compound} @ {conc:.2e} M  -  {out.get('verdict')}", bold=True)
    typer.echo(f"  real effect {_fmt(out.get('real_effect'))} ({out.get('readout')})")
    for m in out.get("modes") or []:
        typer.echo(
            f"  {m['mode']:26} p={m['p_two_sided']:.4f} "
            f"(resolution {m['p_resolution']:.4f})  "
            f"{'beats null' if m['beats_null'] else 'does not beat null'}   z={m['z']:.2f}"
        )
    level = out.get("necessary_information_level") or {}
    if level:
        typer.echo(f"  necessary information level: {level.get('level')} - {level.get('description')}")
    _warnings(out)


@app.command()
def ablation(
    compound: str = typer.Argument("imidacloprid"),
    conc: float = typer.Option(1e-6, "--conc"),
    table: bool = typer.Option(False, "--table", help="Every compound at every concentration."),
    json_out: bool = JSON_OPT,
):
    """What each modelling layer adds: receptor, composition, topology, full."""
    out = _bridge("/api/ablation", {"compound": compound, "conc_M": conc, "table": table})
    if json_out:
        _dump(out)
        return
    _estimate_line(out)
    for level in out.get("level_order") or []:
        effect = (out.get("effects") or {}).get(level)
        typer.echo(f"  {level:22} {_fmt(effect)}")
    for row in (out.get("information_gain") or {}).get("levels") or []:
        typer.echo(
            f"  {row['level']:22} r={_fmt(row.get('pearson_r_vs_full'), 6, 2)} "
            f"rho={_fmt(row.get('spearman_rho_vs_full'), 6, 2)} "
            f"reproduces_full={row.get('reproduces_full')}"
        )
    _warnings(out)


@app.command()
def stability(
    fast: bool = typer.Option(True, "--fast/--full"),
    conc: float = typer.Option(1e-6, "--conc"),
    estimate: bool = typer.Option(False, "--estimate", help="Print the runtime estimate only."),
    json_out: bool = JSON_OPT,
):
    """Every prospective conclusion, re-derived under every admissible rule."""
    out = _bridge(
        "/api/robustness/stability",
        {"fast": fast, "conc_M": conc, "estimate_only": estimate},
    )
    if json_out:
        _dump(out)
        return
    _estimate_line(out)
    if out.get("estimate_only"):
        return
    for row in out.get("rows") or []:
        typer.echo(
            f"  {row['conclusion']:34} retained {row['n_retained']}/{row['n_specs']} "
            f"{'FRAGILE' if row.get('fragile') else 'stable'}"
        )
    _warnings(out)


@app.command()
def uncertainty(
    compound: str = typer.Argument("imidacloprid"),
    conc: float = typer.Option(1e-6, "--conc"),
    n_base: int = typer.Option(32, "--n-base"),
    estimate: bool = typer.Option(False, "--estimate", help="Print the runtime estimate only."),
    json_out: bool = JSON_OPT,
):
    """Global (Sobol) attribution of this model's output variance."""
    out = _bridge(
        "/api/uncertainty/global",
        {"compound": compound, "conc_M": conc, "n_base": n_base, "estimate_only": estimate},
    )
    if json_out:
        _dump(out)
        return
    _estimate_line(out)
    if out.get("estimate_only"):
        return
    typer.echo(f"  output sd {_fmt(out.get('output_sd'))}  n_evaluations {out.get('n_evaluations')}")
    for row in out.get("budget") or []:
        typer.echo(f"  {str(row.get('source')):28} share {_fmt(row.get('share_of_variance'), 8, 3)}")
    _warnings(out)


@app.command()
def voi(
    compound: str = typer.Argument("imidacloprid"),
    conc: float = typer.Option(1e-6, "--conc"),
    n_base: int = typer.Option(32, "--n-base"),
    estimate: bool = typer.Option(False, "--estimate", help="Print the runtime estimate only."),
    json_out: bool = JSON_OPT,
):
    """Which experiment would remove the most model variance."""
    out = _bridge(
        "/api/voi",
        {"compound": compound, "conc_M": conc, "n_base": n_base, "estimate_only": estimate},
    )
    if json_out:
        _dump(out)
        return
    _estimate_line(out)
    if out.get("estimate_only"):
        return
    for row in out.get("rows") or []:
        typer.echo(
            f"  {row.get('rank')}. {str(row.get('factor')):20} "
            f"variance removed {_fmt(row.get('voi_fraction'), 7, 3)}  {row.get('experiment')}"
        )
    _warnings(out)


@app.command()
def claims(
    compound: str = typer.Argument("imidacloprid"),
    conc: float = typer.Option(1e-6, "--conc"),
    markdown: bool = typer.Option(False, "--markdown", help="Render the audit as markdown."),
    json_out: bool = JSON_OPT,
):
    """The dependency chain behind a result: fact, inference and unknown."""
    out = _bridge("/api/claims", {"compound": compound, "conc_M": conc})
    if json_out:
        _dump(out)
        return
    if markdown:
        from flylab.analysis.claims import to_markdown

        typer.echo(to_markdown(out))
        return
    typer.secho(f"claim provenance - {compound} @ {conc:.2e} M", bold=True)
    for link in out.get("chain") or []:
        typer.echo(f"  {link['order']}. {link['step']:26} {link['label']:18} [{link['classification']}]")
        typer.echo(f"     {link['statement']}")
    fiu = out.get("fact_inference_unknown") or {}
    for key, title in (("facts", "FACT"), ("model_inference", "INFERENCE"), ("unknown", "UNKNOWN")):
        typer.echo("")
        typer.secho(title, bold=True)
        for item in fiu.get(key) or []:
            typer.echo(f"  - {item.get('statement')}")
    _warnings(out)


@app.command("reproduce-paper")
def reproduce_paper():
    """Print the command that regenerates every figure and table in the paper."""
    try:
        from scripts.reproduce_paper import main  # type: ignore[import-not-found]
    except ImportError:
        typer.secho(
            "scripts/reproduce_paper.py is not in this checkout yet; it ships with "
            "the paper build. Run it from the repository root once it exists.",
            fg=typer.colors.YELLOW,
        )
        raise typer.Exit(code=2)
    main()


# --------------------------------------------------------------------------
# graph
# --------------------------------------------------------------------------
@graph_app.command("info")
def graph_info(
    graph: str = typer.Option("named", "--graph", help="named | taste_motor | path to JSON"),
    json_out: bool = JSON_OPT,
):
    """Nodes, edges, seeds and transmitter counts of a committed graph."""
    from collections import Counter

    from flylab.circuit.rate import load_graph, resolve_graph

    path = resolve_graph(graph)
    g = load_graph(path)
    nts = Counter((n.get("consensus_nt") or "unclear") for n in g["nodes"])
    payload = {
        "graph": graph,
        "path": str(path),
        "map": g.get("map"),
        "citation": g.get("citation"),
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "hops": g.get("hops"),
        "min_weight": g.get("min_weight"),
        "closure_min_weight": g.get("closure_min_weight"),
        "seeds": {t: len(ids) for t, ids in (g.get("seeds") or {}).items()},
        "transmitters": dict(nts.most_common()),
    }
    if json_out:
        _dump(payload)
        return
    typer.echo(f"{payload['path']}")
    typer.echo(f"map {payload['map']}  hops={payload['hops']}  min_weight={payload['min_weight']}")
    typer.echo(f"{payload['n_nodes']} nodes, {payload['n_edges']} edges")
    typer.echo("seeds: " + ", ".join(f"{t}={n}" for t, n in payload["seeds"].items()))
    typer.echo("transmitters: " + ", ".join(f"{t}={n}" for t, n in payload["transmitters"].items()))


@graph_app.command("impact")
def graph_impact(
    compound: str = "imidacloprid",
    conc: float = 1e-6,
    graph: str = typer.Option("named", "--graph"),
    top: int = 10,
    json_out: bool = JSON_OPT,
):
    """Per-node and per-edge effect of one dose on a graph."""
    from flylab.analysis.impact import summarize_impact

    out = summarize_impact(compound, conc, graph=graph, top_nodes=top, top_edges=top, top_paths=top)
    if json_out:
        _dump(out)
        return
    typer.echo(f"{out['graph']}  {compound} @ {conc:.2e} M")
    typer.echo("  gains " + "  ".join(f"{k}={v:.3f}" for k, v in out["gains"].items()))
    typer.echo(f"{'bodyId':>12} {'type':16} {'nt':16} {'vehicle':>9} {'treated':>9} {'delta':>9}")
    for row in out["top_nodes"][:top]:
        typer.echo(
            f"{row.get('bodyId'):>12} {str(row.get('type'))[:16]:16} {str(row.get('nt'))[:16]:16} "
            f"{_fmt(row.get('rate_vehicle'), 9)} {_fmt(row.get('rate_treated'), 9)} {_fmt(row.get('delta_hz'), 9)}"
        )
    _warnings(out)


# --------------------------------------------------------------------------
# experiments
# --------------------------------------------------------------------------
EXAMPLE_DESIGN = """\
# FlyLab experiment design (feed to: flylab experiment run DESIGN.yaml)
assay: subgraph            # subgraph | spiking | taste | taste_map
compounds:
  - imidacloprid
  - nicotine
  - fipronil
concs_M: [1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
replicates: 3
seed: 0
readouts: [mn9_hz, dnp01_hz, mean_hz, g_ach, g_gaba]
jitter_log10: 0.3          # 0 disables the library Monte-Carlo
include_vehicle: true
graph: named               # named | taste_motor
"""


@experiment_app.command("example")
def experiment_example():
    """Print a starter design YAML you can redirect into a file."""
    typer.echo(EXAMPLE_DESIGN, nl=False)


@experiment_app.command("run")
def experiment_run(
    design: Path = typer.Argument(..., help="Design YAML or JSON file (a spec file also works)."),
    out: Path | None = typer.Option(None, "--out", help="Write the results CSV here."),
    json_out_path: Path | None = typer.Option(None, "--json", help="Write the full JSON result here."),
    outdir: Path | None = typer.Option(
        None, "--outdir", help="Write a full run directory instead; identical to `flylab run`."
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate and estimate only; run nothing."),
):
    """Run a batch design and print (or write) the results table.

    A thin alias of `flylab run`: with --outdir or --dry-run it is that command
    exactly, and a v0.5 design file is a valid spec (``assay``, ``concs_M`` and
    ``compound`` are accepted aliases). Without them it keeps the v0.5
    behaviour of emitting the flat results table.
    """
    if outdir is not None or dry_run:
        run_spec_command(spec=design, outdir=outdir, dry_run=dry_run, json_out=False)
        return

    from flylab.assays.experiment import design_from_yaml, run_experiment

    spec = design_from_yaml(design)
    result = run_experiment(spec)
    if out:
        Path(out).write_text(result["csv"])
        typer.echo(f"csv  -> {out}  ({result['n_rows']} drug rows)")
    if json_out_path:
        Path(json_out_path).write_text(json.dumps(result, indent=2, default=str))
        typer.echo(f"json -> {json_out_path}")
    if not out and not json_out_path:
        typer.echo(result["csv"], nl=False)
    typer.echo(
        f"# {result['n_rows']} rows + {len(result['vehicle_rows'])} vehicle rows, "
        f"{len(result['summary'])} groups",
        err=True,
    )
    for w in result.get("warnings") or []:
        typer.secho(f"# ! {w}", fg=typer.colors.YELLOW, err=True)


# --------------------------------------------------------------------------
# specs: one YAML file describes a whole run
# --------------------------------------------------------------------------
@app.command("run")
def run_spec_command(
    spec: Path = typer.Argument(..., help="Experiment spec YAML (see `flylab spec-schema --example`)."),
    outdir: Path | None = typer.Option(
        None, "--outdir", help="Run directory to write. Defaults to runs/<spec name>."
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Validate, print the resolved spec and a runtime estimate; run nothing."
    ),
    json_out: bool = JSON_OPT,
):
    """Run a declarative experiment spec into a self-describing run directory."""
    from flylab.spec import SpecError, load_spec, resolved_spec_dict, run_spec

    try:
        parsed = load_spec(spec)
    except SpecError as exc:
        typer.secho(f"{spec}: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2)

    target = Path(outdir) if outdir else Path("runs") / parsed.name
    try:
        out = run_spec(
            parsed,
            target,
            dry_run=dry_run,
            progress=None if json_out else (lambda msg: typer.echo(msg)),
        )
    except SpecError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2)

    if json_out:
        _dump(out)
        return

    if dry_run:
        import yaml

        typer.secho(f"resolved spec ({spec})", bold=True)
        typer.echo(yaml.safe_dump(json.loads(json.dumps(resolved_spec_dict(parsed), default=str)), sort_keys=False))
        est = out["runtime_estimate"]
        typer.secho(
            f"estimate: ~{est['estimate_human']} "
            f"({est['n_assay_runs']} assay runs"
            + (f" + {', '.join(sorted(est['analyses']))}" if est["analyses"] else "")
            + ")",
            fg=typer.colors.BLUE,
        )
        for name, detail in sorted(est["analyses"].items()):
            typer.echo(f"  {name:12} ~{detail['seconds']} s")
        typer.echo(f"would write -> {target}/  ({', '.join(out['would_write'])})")
        typer.echo(f"spec_sha256 {out['spec_sha256']}")
        for w in out.get("warnings") or []:
            typer.secho(f"  ! {w}", fg=typer.colors.YELLOW)
        return

    typer.secho(f"run -> {target}", bold=True)
    typer.echo(f"  rows      {out['n_rows']}")
    typer.echo(f"  notebooks {len(out['notebooks'])}")
    typer.echo(f"  analyses  {', '.join(out['analyses_completed']) or 'none'}")
    typer.echo(f"  cards     {len([c for c in out['cards'] if c.endswith('.json')])}")
    typer.echo(f"  digest    {out['determinism']['digest']}")
    typer.echo(f"  total     {out['timings_s'].get('total_s')} s")
    _warnings(out)


@app.command("spec-schema")
def spec_schema_command(
    example: bool = typer.Option(False, "--example", help="Print a starter spec YAML instead."),
    json_out: bool = JSON_OPT,
):
    """The machine-readable schema for an experiment spec."""
    from flylab.spec import EXAMPLE_SPEC, spec_schema

    if example:
        typer.echo(EXAMPLE_SPEC, nl=False)
        return
    schema = spec_schema()
    if json_out:
        _dump(schema)
        return
    typer.secho(f"{schema['title']} v{schema['flylab_spec_version']}", bold=True)
    for name, prop in schema["properties"].items():
        default = prop.get("default")
        typer.echo(f"  {name:22} {str(prop.get('type')):8} default={json.dumps(default, default=str)}")
        typer.echo(f"      {prop.get('description')}")
        if prop.get("enum"):
            typer.echo(f"      one of: {', '.join(str(v) for v in prop['enum'])}")
    typer.secho("analyses", bold=True)
    for name, meta in schema["x-analyses"].items():
        typer.echo(f"  {name:12} -> {meta['artifact']}  ({meta['module']})")
        for opt, text in meta["options"].items():
            typer.echo(f"      {opt:16} {text}")


@app.command("card")
def card_command(
    compound: str = typer.Argument("imidacloprid", help="Compound, or -- with --run -- ignored."),
    conc: float = typer.Option(1e-6, "--conc", help="Free concentration in molar."),
    engine: str = typer.Option("rate", "--engine", help="rate | lif"),
    graph: str = typer.Option("named", "--graph"),
    run: Path | None = typer.Option(None, "--run", help="Read the cards of an existing run directory."),
    json_out: bool = typer.Option(False, "--json", help="Print the card(s) as JSON."),
    markdown: bool = typer.Option(False, "--markdown", help="Print the card(s) as Markdown."),
):
    """A claim card: what is claimed, what backs it, and what was assumed."""
    from flylab.report.card import cards_for_run, claim_card, to_markdown

    if run is not None:
        try:
            cards = cards_for_run(run)
        except FileNotFoundError as exc:
            typer.secho(str(exc), fg=typer.colors.RED, err=True)
            raise typer.Exit(code=2)
        if not cards:
            typer.secho(f"{run}: no cards and no notebooks to build them from", fg=typer.colors.YELLOW)
            raise typer.Exit(code=2)
    else:
        from flylab.assays.ensemble import run_assay

        assay = "spiking" if engine in ("lif", "spiking") else "subgraph"
        notebook = run_assay(assay, compound, conc, graph=graph)
        vehicle = run_assay(assay, None, 0.0, graph=graph)
        cards = [
            claim_card(
                notebook,
                compound=compound,
                conc_M=conc,
                engine="lif" if assay == "spiking" else "rate",
                graph=graph,
                vehicle=vehicle,
            )
        ]

    if json_out:
        _dump(cards[0] if len(cards) == 1 else cards)
        return
    if markdown:
        typer.echo("\n".join(to_markdown(c) for c in cards), nl=False)
        return
    for card in cards:
        subject = card.get("subject") or {}
        typer.secho(
            f"claim card - {subject.get('compound')} @ "
            f"{(subject.get('concentration_M') or 0):.2e} M ({subject.get('engine')})",
            bold=True,
        )
        typer.echo(f"  {card['claim']['sentence']}")
        for key in card.get("sections") or []:
            section = card.get(key) or {}
            mark = " " if section.get("assessed", True) else "!"
            text = section.get("statement") or section.get("sentence") or ""
            typer.echo(f"  {mark} {key:36} {str(text)[:96]}")
    _warnings(cards[0])


@app.command("manifest")
def manifest_command(
    write: bool = typer.Option(False, "--write", help="Write the repository artifact manifest."),
    check: bool = typer.Option(False, "--check", help="Verify the committed manifest; exit 2 on mismatch."),
    lock: bool = typer.Option(False, "--lock", help="Also write requirements-lock.txt."),
    path: Path | None = typer.Option(None, "--path", help="Manifest path (default artifact-manifest.json)."),
    json_out: bool = JSON_OPT,
):
    """Compute or verify the repository artifact manifest."""
    manifest_mod = _manifest_module()
    target = Path(path) if path else Path(manifest_mod.MANIFEST_PATH)

    if lock:
        written = manifest_mod.write_lockfile()
        typer.echo(f"lock -> {written['path']}  ({written['n_packages']} packages)")

    if write:
        payload = manifest_mod.build_manifest()
        manifest_mod.write_manifest(payload, target)
        if json_out:
            _dump(payload)
        else:
            typer.secho(f"manifest -> {target}", bold=True)
            for section, entry in sorted(payload["sections"].items()):
                typer.echo(
                    f"  {section:13} {entry['n_entries']:4d} entries  {entry['sha256'][:12]}"
                    + ("" if entry.get("fatal") else "  (informational)")
                )
        if not check:
            return

    if check or not write:
        if not target.exists():
            typer.secho(
                f"{target} does not exist yet. Create it with `flylab manifest --write`.",
                fg=typer.colors.RED,
                err=True,
            )
            raise typer.Exit(code=2)
        diff = manifest_mod.verify_manifest(target)
        if json_out:
            _dump(diff)
        else:
            typer.secho(f"manifest check - {target}", bold=True)
            for section, result in diff["sections"].items():
                status = "ok" if result["ok"] else "MISMATCH"
                colour = typer.colors.GREEN if result["ok"] else typer.colors.RED
                typer.secho(f"  {section:12} {status:9} {result['summary']}", fg=colour)
            for note in diff.get("notes") or []:
                typer.secho(f"  note: {note}", fg=typer.colors.YELLOW)
        if not diff["ok"]:
            raise typer.Exit(code=2)


def _manifest_module():
    """Import ``scripts/manifest.py`` from an installed package or a checkout."""
    import importlib
    import importlib.util
    import sys

    try:
        return importlib.import_module("scripts.manifest")
    except Exception:
        pass
    here = Path(__file__).resolve().parents[1] / "scripts" / "manifest.py"
    if not here.exists():
        typer.secho(
            "scripts/manifest.py is not in this checkout; the artifact manifest is a "
            "repository tool and is not shipped in the wheel.",
            fg=typer.colors.YELLOW,
            err=True,
        )
        raise typer.Exit(code=2)
    spec_obj = importlib.util.spec_from_file_location("flylab_scripts_manifest", here)
    module = importlib.util.module_from_spec(spec_obj)  # type: ignore[arg-type]
    sys.modules["flylab_scripts_manifest"] = module
    spec_obj.loader.exec_module(module)  # type: ignore[union-attr]
    return module


# --------------------------------------------------------------------------
# data / serving
# --------------------------------------------------------------------------
@app.command("download-malecns")
def download_malecns(full: bool = typer.Option(False, "--full")):
    """Download the MaleCNS atlas (and with --full, the 1.1 GB weight matrix)."""
    typer.echo(str(download_malecns_files(kind="full" if full else "atlas")))


@app.command("extract-subgraph")
def extract_subgraph(
    hops: int = 1,
    min_weight: int = 5,
    types: str = typer.Option("MN9,DNp01", help="Comma-separated MaleCNS type names"),
    closure_min_weight: int | None = typer.Option(
        None,
        "--closure-min-weight",
        help="Weight floor for the induced edges between neighborhood members "
        "(defaults to --min-weight). Raise it to keep the dense taste_motor JSON small.",
    ),
    out: Path = typer.Option(
        Path("data/derived/malecns_named_neighborhood.json"),
        "--out",
        help="Destination JSON path.",
    ),
):
    """Cut a named-cell neighborhood out of the MaleCNS weight matrix."""
    from flylab.maps.extract import extract_neighborhood

    type_list = tuple(t.strip() for t in types.split(",") if t.strip())
    payload = extract_neighborhood(
        types=type_list,
        hops=hops,
        min_weight=min_weight,
        closure_min_weight=closure_min_weight,
        out=Path(out),
    )
    typer.echo(f"seeds={payload['seeds']}")
    typer.echo(f"nodes={payload['n_nodes']} edges={payload['n_edges']} -> {payload['path']}")


# --------------------------------------------------------------------------
# tutorial
# --------------------------------------------------------------------------
# Every lesson here is the executable half of docs/TUTORIAL.md.  A step is a
# real invocation of a real command: the runner drives this same Typer app, so
# nothing in a lesson can drift away from what `flylab <command>` does.  The
# `not_` line of each lesson is the half a newcomer skips, which is why it is
# printed after the output rather than before it.
#
# Nothing here reaches the network.  `n` values are deliberately below the
# paper's; each lesson names the paper-scale command in its notes.
TUTORIAL_DOC = "docs/TUTORIAL.md"
INTERPRETATION_DOC = "docs/INTERPRETATION.md"


def _tutorial_own_compound() -> None:
    """Lesson 7's step: add a row, then try to type it wrongly.

    Kept in Python because the library is an input, not a command-line flag:
    ``compare_compound(..., library=...)`` is the supported way to score a
    compound the shipped YAML does not contain.
    """
    from flylab.pharm import evidence as ev
    from flylab.pharm.occupancy import compare_compound, load_library

    lib = load_library()
    lib["compounds"]["my_neonic"] = {
        "name": "My neonicotinoid",
        "class": "neonicotinoid",
        "receptors": {
            "insect_nAChR": {
                "param_type": "EC50",
                "value_M": 4.0e-8,
                "n": 1.1,
                "direction": "agonist",
                "relation": "exact_compound_exact_receptor_exact_species",
                "species": "Drosophila melanogaster",
                "source": "my lab, two-electrode voltage clamp on Dalpha1/Dbeta1",
                "evidence_tier": "literature_order",
            },
            "vertebrate_nAChR_a4b2": {
                "param_type": "unknown",
                "value_M": None,
                "direction": "none",
                "relation": "unsupported",
                "source": "no source measured this compound at this receptor",
                "evidence_tier": "class_placeholder",
            },
        },
    }
    result = compare_compound("my_neonic", 1e-6, library=lib)
    typer.echo(f"{'receptor':24} {'engagement':>10}  {'model':28} distance")
    for row in result["receptors"]:
        eng = row.get("engagement")
        typer.echo(
            f"{row['receptor']:24} "
            f"{('       n/a' if eng is None else f'{eng:10.3f}')}  "
            f"{str(row.get('engagement_model')):28} {row.get('evidence_distance')}"
        )

    typer.echo("\nnow type the same row wrongly:")
    for pt, model, relation, why in (
        ("EC50", "binding_occupancy", "exact_compound_exact_receptor_exact_species",
         "an EC50 relabelled as an occupancy"),
        ("Kd", "binding_occupancy", "exact_compound_exact_receptor_other_species",
         "a binding constant measured in another species"),
        ("unknown", "functional_engagement", "unsupported",
         "a placeholder asked for a number"),
    ):
        typer.echo(f"\n  {why}:")
        try:
            ev.check_transformation(pt, model, relation)
        except ev.EvidenceTypeError as exc:
            typer.secho(f"    EvidenceTypeError: {exc}", fg=typer.colors.RED)
        else:  # pragma: no cover - would be a regression in evidence.py
            typer.secho("    ACCEPTED - this is a bug in flylab/pharm/evidence.py", fg=typer.colors.RED)

    typer.echo("\n  the vocabulary is closed too:")
    for call, arg in (("as_param_type", "Ki50"), ("as_relation", "measured_in_our_lab")):
        try:
            getattr(ev, call)(arg)
        except ev.EvidenceTypeError as exc:
            typer.secho(f"    EvidenceTypeError: {exc}", fg=typer.colors.RED)


#: The lessons, in the order docs/TUTORIAL.md teaches them.  ``steps`` are run
#: end to end by ``flylab tutorial <n>``; ``argv`` is fed to this same app.
TUTORIAL_LESSONS: tuple[dict[str, Any], ...] = (
    {
        "id": 1,
        "key": "dashboard",
        "title": "The dashboard for one compound at one dose",
        "teaches": "what each headline tile is, and which of them is a measurement",
        "not_": (
            "The four tiles are three model outputs and one evidence census. None of "
            "them is an observation, and the engagement tile is not an occupancy."
        ),
        "steps": [
            {
                "argv": ["dashboard", "imidacloprid", "--conc", "1e-6"],
                "note": (
                    "Insect engagement reads 1.000 at insect_nAChR_beta1 - and the line "
                    "beside it says binding_engagement_proxy, E1 (cross-species). A full "
                    "bar on an aphid binding constant is still an extrapolation onto "
                    "Drosophila. See docs/INTERPRETATION.md, 'Engagement is not occupancy'."
                ),
            },
            {
                "argv": ["occupancy", "imidacloprid", "--conc", "1e-6"],
                "note": (
                    "Two of the six rows print n/a and 'not modelled'. They are excluded "
                    "from every number above, and they are not zeros."
                ),
            },
        ],
    },
    {
        "id": 2,
        "key": "local",
        "title": "The same numbers from the CLI, the server and the browser",
        "teaches": "that the three front ends call one implementation",
        "not_": (
            "Identical numbers prove the transports agree. They say nothing about "
            "whether the model is right."
        ),
        "steps": [
            {
                "argv": ["assay-subgraph", "--compound", "imidacloprid", "--conc", "1e-6"],
                "note": "mean_hz here is the number the dashboard's circuit card prints.",
            },
            {
                "argv": ["meta"],
                "note": (
                    "The library sha256 and the map id are what make two runs comparable. "
                    "Quote them with any number you take out of the bench."
                ),
            },
        ],
    },
    {
        "id": 3,
        "key": "dose",
        "title": "A dose-response, and the parameter that stops mattering",
        "teaches": "why a saturating dose makes the cited potency irrelevant",
        "not_": (
            "The EC50 is irrelevant *at this dose*. Between 10 nM and 100 nM it is the "
            "steepest thing in the model. 'Insensitive to potency' is a statement about "
            "one point on the curve, never about the model."
        ),
        "steps": [
            {"argv": ["assay-subgraph", "--compound", "imidacloprid", "--conc", "1e-8"]},
            {"argv": ["assay-subgraph", "--compound", "imidacloprid", "--conc", "1e-7"]},
            {"argv": ["assay-subgraph", "--compound", "imidacloprid", "--conc", "1e-6"]},
            {
                "argv": ["assay-subgraph", "--compound", "imidacloprid", "--conc", "1e-5"],
                "note": (
                    "1e-6 and 1e-5 agree to every printed digit: the receptor is full and "
                    "the gain patch has hit its floor."
                ),
            },
            {
                "argv": ["sensitivity", "--compound", "imidacloprid", "--conc", "1e-6"],
                "note": (
                    "A two-fold error in the EC50 moves mean_hz by 0.000 Hz; the asserted "
                    "gain coefficient moves it by 2.24 Hz. That ordering is why the "
                    "uncertainty budget puts potency at the estimator's noise floor and "
                    "the gain rule at the top - `flylab voi imidacloprid --conc 1e-6`."
                ),
            },
        ],
    },
    {
        "id": 4,
        "key": "compare",
        "title": "Two compounds whose selectivity orderings disagree",
        "teaches": "that a receptor ratio does not predict a circuit window",
        "not_": (
            "Neither index is a safety margin. A circuit selectivity index says only "
            "that this simulation moves before the teaching library's vertebrate "
            "receptor fills."
        ),
        "steps": [
            {
                "argv": ["compare", "imidacloprid", "deltamethrin", "--conc", "1e-6",
                         "--no-occupancy", "--n", "20"],
                "note": (
                    "Imidacloprid wins on receptor selectivity (2.70 log10 against 2.07) "
                    "and loses on circuit selectivity (1.79 against 2.28). The two "
                    "quantities are computed from different things and rank differently."
                ),
            },
        ],
    },
    {
        "id": 5,
        "key": "dependence",
        "title": "Did the connectome matter? - and the answer changes with the cut",
        "teaches": "the three verdicts, and that a verdict belongs to a substrate",
        "not_": (
            "A large permutation probability licenses 'not distinguishable from this "
            "null ensemble at n shuffles' and nothing else. Equivalence needs the gap "
            "from the null median to fall below the prespecified margin, and the "
            "result says when it does. And neither verdict below is the truth about "
            "the compound: the answer moves with the cut AND with its size, so the "
            "real question is where it settles - which is still being measured. "
            "Never quote a verdict without the cut, its size and its mean degree."
        ),
        "steps": [
            {
                "argv": ["dependence", "imidacloprid", "--conc", "1e-6", "--n", "200",
                         "--graph", "named"],
                "note": (
                    "composition-dominated on the `named` cut. About 85% of that cut's "
                    "edges land on four seed cells, so a degree-preserving rewire is "
                    "nearly the identity: a negative result here is weak evidence."
                ),
            },
            {
                "argv": ["dependence", "imidacloprid", "--conc", "1e-6", "--n", "40",
                         "--graph", "taste_motor"],
                "note": (
                    "Same compound, same dose, same readout, opposite verdict. The "
                    "paper-scale command is the same line with --n 1000. The scale "
                    "ladder (flylab.analysis.scale.available_cuts) moves it again with "
                    "size alone: composition-dominated at ~1100 cells, distinguishable "
                    "from every null at ~5000. See docs/INTERPRETATION.md section 4."
                ),
            },
        ],
    },
    {
        "id": 6,
        "key": "spec",
        "title": "A batch run from a spec file, and the claim card it writes",
        "teaches": "how a result arrives with everything needed to disbelieve it",
        "not_": (
            "Sections 6 and 7 of a card say 'not assessed in this run' when the spec "
            "did not ask for them. That is an unanswered question, not an absent "
            "concern."
        ),
        "steps": [
            {
                "argv": ["spec-schema", "--example"],
                "note": "Write that to a file and `flylab run it.yaml --dry-run` validates it.",
            },
            {
                "argv": ["card", "fipronil", "--conc", "1e-6", "--markdown"],
                "note": (
                    "The same nine sections a run directory writes to cards/. Section 3 "
                    "is where the evidence typing shows: fipronil's insect_RDL row is an "
                    "IC50 measured in Drosophila (functional engagement, E0), its "
                    "vertebrate GABA-A row an IC50 in human recombinant receptors "
                    "(functional engagement proxy, E1)."
                ),
            },
        ],
    },
    {
        "id": 7,
        "key": "library",
        "title": "Bringing your own compound, and being refused",
        "teaches": "the parameter type, the evidence distance and the source string",
        "not_": (
            "The type system checks that you described your measurement honestly. It "
            "cannot check that the measurement is right, and it never invents a value "
            "for a receptor you left out."
        ),
        "steps": [
            {
                "python": _tutorial_own_compound,
                "label": "python - add a row, then mistype it",
                "note": (
                    "A row is admitted on (param_type, relation). Nothing you write in "
                    "`source` can raise the engagement model, and nothing can turn an "
                    "unsupported row into a small number."
                ),
            },
        ],
    },
)


def _tutorial_lesson(selector: str) -> dict[str, Any]:
    key = str(selector).strip().lower()
    for lesson in TUTORIAL_LESSONS:
        if key == str(lesson["id"]) or key == lesson["key"]:
            return lesson
    valid = ", ".join(f"{le['id']} ({le['key']})" for le in TUTORIAL_LESSONS)
    typer.secho(f"no tutorial lesson {selector!r}. Lessons: {valid}", fg=typer.colors.RED, err=True)
    raise typer.Exit(code=2)


def _tutorial_command_line(step: dict[str, Any]) -> str:
    if "argv" in step:
        return "flylab " + " ".join(step["argv"])
    return str(step.get("label") or "python")


def _tutorial_run_step(step: dict[str, Any]) -> None:
    typer.secho("$ " + _tutorial_command_line(step), bold=True)
    if "python" in step:
        step["python"]()
        return
    from typer.main import get_command

    try:
        get_command(app).main(
            args=list(step["argv"]), prog_name="flylab", standalone_mode=False
        )
    except SystemExit:  # pragma: no cover - click raises typer.Exit instead
        pass
    except click_exceptions.Exit as exc:  # a command that ended early, not a crash
        if int(getattr(exc, "exit_code", 0) or 0) not in (0, 2):
            raise


def _tutorial_run(lesson: dict[str, Any]) -> None:
    typer.secho(f"\nLesson {lesson['id']} - {lesson['title']}", bold=True)
    typer.secho(f"  teaches: {lesson['teaches']}", fg=typer.colors.CYAN)
    for step in lesson["steps"]:
        typer.echo("")
        _tutorial_run_step(step)
        if step.get("note"):
            typer.secho("  -> " + step["note"], fg=typer.colors.CYAN)
    typer.echo("")
    typer.secho("  what this does NOT tell you", bold=True)
    typer.secho("  " + lesson["not_"], fg=typer.colors.YELLOW)
    typer.echo(
        f"\n  full lesson: {TUTORIAL_DOC}   how to read the output: {INTERPRETATION_DOC}"
    )


@app.command("tutorial")
def tutorial_command(
    lesson: str = typer.Argument(
        None, help="Lesson number or key. Omit it (or pass --list) to see them all."
    ),
    list_lessons: bool = typer.Option(False, "--list", help="Print the lesson list and stop."),
    json_out: bool = JSON_OPT,
):
    """Guided lessons: print the list, or run one end to end with its caveats.

    ``flylab tutorial --list`` is the index of docs/TUTORIAL.md;
    ``flylab tutorial 3`` runs lesson 3's real commands and prints the real
    output with the interpretation note each number needs.  Nothing here
    touches the network.
    """
    if json_out and not lesson:
        _dump(
            [
                {
                    "id": le["id"],
                    "key": le["key"],
                    "title": le["title"],
                    "teaches": le["teaches"],
                    "does_not_tell_you": le["not_"],
                    "commands": [_tutorial_command_line(s) for s in le["steps"]],
                }
                for le in TUTORIAL_LESSONS
            ]
        )
        return
    if list_lessons or not lesson:
        typer.secho("FlyLab tutorial - " + TUTORIAL_DOC, bold=True)
        for le in TUTORIAL_LESSONS:
            typer.echo(f"  {le['id']}. {le['title']}  [{le['key']}]")
            typer.secho(f"     {le['teaches']}", fg=typer.colors.CYAN)
        typer.echo("\n  run one:   flylab tutorial 3")
        typer.echo(f"  read it:   {TUTORIAL_DOC}")
        typer.echo(f"  read the outputs: {INTERPRETATION_DOC}")
        return
    chosen = _tutorial_lesson(lesson)
    if json_out:
        _dump(
            {
                "id": chosen["id"],
                "key": chosen["key"],
                "title": chosen["title"],
                "teaches": chosen["teaches"],
                "does_not_tell_you": chosen["not_"],
                "commands": [_tutorial_command_line(s) for s in chosen["steps"]],
            }
        )
        return
    _tutorial_run(chosen)


@app.command()
def serve(
    host: str = "127.0.0.1",
    port: int = 8765,
    open_browser: bool = typer.Option(False, "--open", help="Open the bench in a browser."),
):
    """Serve the bench UI and HTTP API."""
    import uvicorn

    url = f"http://{host}:{port}"
    typer.echo(f"FlyLab bench -> {url}")
    if open_browser:
        import threading
        import webbrowser

        threading.Timer(1.2, lambda: webbrowser.open(url)).start()
    uvicorn.run("flylab.server:app", host=host, port=port, reload=False)


if __name__ == "__main__":
    app()
