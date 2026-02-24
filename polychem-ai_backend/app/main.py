from contextlib import asynccontextmanager
import os

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse

from app.schemas import RecommendRequest, RecommendResponse
from app.store import load_dataset
from app.core import recommend_new_compound
from app.settings import STATIC_DIR, COMPOUNDS_DIR
from app.cache import get_history


@asynccontextmanager
async def lifespan(app: FastAPI):
    print("Loading dataset...")
    load_dataset()
    print("Dataset loaded!")
    yield


app = FastAPI(
    title="Novachem AI Backend",
    version="1.0.0",
    lifespan=lifespan,
)

# =========================================================
# Static files (Koyeb-safe: /tmp writable)
# =========================================================
os.makedirs(STATIC_DIR, exist_ok=True)
os.makedirs(COMPOUNDS_DIR, exist_ok=True)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# =========================================================
# CORS
# - Untuk produksi: isi FRONTEND_URL di environment (lebih aman)
# - Kalau belum ada, fallback ke localhost dev
# =========================================================
frontend_url = os.getenv("FRONTEND_URL")  # contoh: https://your-frontend.koyeb.app

allow_origins = []
if frontend_url:
    allow_origins.append(frontend_url)
    print(f"✅ CORS: Added frontend_url from env: {frontend_url}")
else:
    print(f"⚠️ WARNING: FRONTEND_URL not set. Using fallback origins.")

# Jika tidak ada origins, allow localhost untuk dev
if not allow_origins:
    allow_origins = [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "http://localhost:5174",
        "http://127.0.0.1:5174",
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost",
        "http://127.0.0.1",
        "*",  # fallback untuk debug
    ]

print(f"CORS allow_origins: {allow_origins}")

# Configure CORS properly - MUST be added after app creation but before endpoints
app.add_middleware(
    CORSMiddleware,
    allow_origins=allow_origins,
    allow_credentials=False,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
    expose_headers=["*"],
    max_age=3600,
)

# =========================================================
# API Endpoints
# =========================================================

# Handle CORS preflight requests explicitly
@app.options("/{path_name:path}")
def handle_preflight(path_name: str):
    """Handle CORS preflight requests"""
    return {"message": "OK"}


@app.post("/predict", response_model=RecommendResponse)
def predict(req: RecommendRequest):
    try:
        result = recommend_new_compound(req.smiles)
    except Exception as e:
        print(f"Error di pipeline: {repr(e)}")
        raise HTTPException(status_code=503, detail="Pipeline unavailable")

    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])

    return result


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/history")
def history():
    data = get_history()

    # Backfill defaults biar history aman walau schema berubah
    for item in data:
        res = item.get("result", {})

        nc = res.get("new_compound", {})
        nc.setdefault("formula", "")
        nc.setdefault("molecular_weight", 0.0)
        nc.setdefault("tg", 0.0)
        nc.setdefault("tg_justification", "")
        nc.setdefault("pid", "")
        nc.setdefault("polymer_class", "")
        res["new_compound"] = nc

        scs = res.get("similar_compounds", [])
        for sc in scs:
            sc.setdefault("name", "")
            sc.setdefault("formula", "")
            sc.setdefault("molecular_weight", 0.0)
            sc.setdefault("tg", 0.0)
            sc.setdefault("pid", "")
            sc.setdefault("polymer_class", "")

            score = sc.get("similarity_score", 0.0)
            sc.setdefault("similarity_percent", float(score) * 100.0)

        res["similar_compounds"] = scs
        item["result"] = res

    return JSONResponse(content=jsonable_encoder(data))