from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from flylab.assays.taste import dose_response, library_keys, run_taste_assay
from flylab.pharm.occupancy import compare_compound

STATIC = Path(__file__).parent / "static"

app = FastAPI(title="FlyLab", version="0.2.0")
app.mount("/static", StaticFiles(directory=STATIC), name="static")


class AssayRequest(BaseModel):
    compound: str | None = None
    conc_M: float = Field(0.0, ge=0)
    sugar_hz: float = Field(150.0, ge=0)
    bitter_hz: float = Field(0.0, ge=0)


@app.get("/")
def index() -> FileResponse:
    return FileResponse(STATIC / "index.html")


@app.get("/api/drugs")
def drugs() -> dict:
    return {"compounds": library_keys()}


@app.get("/api/occupancy")
def occupancy(compound: str, conc_M: float) -> dict:
    return compare_compound(compound, conc_M)


@app.post("/api/assay/taste")
def assay_taste(req: AssayRequest) -> dict:
    return run_taste_assay(
        compound=req.compound,
        conc_M=req.conc_M,
        sugar_hz=req.sugar_hz,
        bitter_hz=req.bitter_hz,
    )


@app.get("/api/assay/taste/dose-response")
def assay_dose(compound: str) -> dict:
    concs = [10**e for e in range(-10, -4)]
    return {"compound": compound, "points": dose_response(compound, concs)}
