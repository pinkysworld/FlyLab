from __future__ import annotations
from pathlib import Path
from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from flylab.assays.taste import dose_response, library_keys, run_taste_assay
from flylab.assays.wholens import run_wholens_assay
from flylab.pharm.occupancy import compare_compound

STATIC = Path(__file__).parent / "static"
app = FastAPI(title="FlyLab", version="0.3.0")
app.mount("/static", StaticFiles(directory=STATIC), name="static")

class AssayRequest(BaseModel):
    compound: str | None = None
    conc_M: float = Field(0.0, ge=0)
    sugar_hz: float = Field(150.0, ge=0)
    bitter_hz: float = Field(0.0, ge=0)

@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")

@app.get("/api/drugs")
def drugs():
    return {"compounds": library_keys()}

@app.get("/api/occupancy")
def occupancy(compound: str, conc_M: float):
    return compare_compound(compound, conc_M)

@app.post("/api/assay/taste")
def assay_taste(req: AssayRequest):
    return run_taste_assay(req.compound, req.conc_M, req.sugar_hz, req.bitter_hz)

@app.get("/api/assay/taste/dose-response")
def assay_dose(compound: str):
    return {"compound": compound, "points": dose_response(compound, [10**e for e in range(-10, -4)])}

@app.post("/api/assay/cns")
def assay_cns(req: AssayRequest):
    return run_wholens_assay(req.compound, req.conc_M)
