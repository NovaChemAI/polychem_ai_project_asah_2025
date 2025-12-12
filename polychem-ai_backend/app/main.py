from contextlib import asynccontextmanager
import os
from typing import List

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse

from app.schemas import RecommendRequest, RecommendResponse, HistoryItemOut
from app.store import load_dataset
from app.core import recommend_new_compound
from app.settings import STATIC_DIR
from app.cache import get_history


@asynccontextmanager
async def lifespan(app: FastAPI):
    load_dataset()
    yield


app = FastAPI(
    title="Novachem AI Backend",
    version="1.0.0",
    lifespan=lifespan,
)

os.makedirs(STATIC_DIR, exist_ok=True)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.post("/recommend", response_model=RecommendResponse)
async def recommend(req: RecommendRequest):
    try:
        result = recommend_new_compound(req.smiles)
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Pipeline error: {str(e)}")

    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])

    return result


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/history", response_model=None)
def history():
    data = get_history()

    # backfill
    for item in data:
        res = item.get("result", {})

        # new_compound
        nc = res.get("new_compound", {})
        nc.setdefault("formula", "")
        nc.setdefault("molecular_weight", 0.0)
        nc.setdefault("tg", 0.0)
        nc.setdefault("pid", "")
        nc.setdefault("polymer_class", "")
        res["new_compound"] = nc

        # similar_compounds
        scs = res.get("similar_compounds", [])
        for sc in scs:
            sc.setdefault("name", "")
            sc.setdefault("formula", "")
            sc.setdefault("molecular_weight", 0.0)
            sc.setdefault("tg", 0.0)
            sc.setdefault("pid", "")
            sc.setdefault("polymer_class", "")
            sc.setdefault(
                "similarity_percent",
                float(sc.get("similarity_score", 0.0)) * 100
            )
        res["similar_compounds"] = scs

        item["result"] = res

    return JSONResponse(content=jsonable_encoder(data))