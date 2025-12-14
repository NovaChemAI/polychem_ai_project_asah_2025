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

# Pastikan import ini sesuai dengan struktur folder Anda
from app.schemas import RecommendRequest, RecommendResponse
from app.store import load_dataset
from app.core import recommend_new_compound
from app.settings import STATIC_DIR
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

# Setup Static Files
os.makedirs(STATIC_DIR, exist_ok=True)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# --- PERBAIKAN CORS ---
# Hapus definisi ganda. Masukkan semua port frontend di sini.
origins = [
    "http://localhost:5173",    # Default Vite
    "http://127.0.0.1:5173",    # Default Vite IP
    "http://localhost:3000",    # React CRA / Next
    "http://127.0.0.1:3000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],      # Gunakan list spesifik, jangan "*" jika pakai credentials
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
# ----------------------


# --- PERBAIKAN ENDPOINT ---
# 1. Ganti '/recommend' jadi '/predict' agar cocok dengan frontend Anda
# 2. Hapus 'async' agar server tidak nge-freeze saat AI berpikir
@app.post("/predict", response_model=RecommendResponse)
def predict(req: RecommendRequest):
    try:
        # Proses AI berat berjalan di sini
        result = recommend_new_compound(req.smiles)
    except Exception as e:
        print(f"Error di pipeline: {e}") # Print error ke terminal backend
        raise HTTPException(status_code=503, detail=f"Pipeline error: {str(e)}")

    if "error" in result:
        # Ini akan mengirim pesan error seperti "SMILES tidak valid" ke frontend
        raise HTTPException(status_code=400, detail=result["error"])

    return result


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/history")
def history():
    data = get_history()

    # Backfill logic
    for item in data:
        res = item.get("result", {})

        # new_compound defaults
        nc = res.get("new_compound", {})
        nc.setdefault("formula", "")
        nc.setdefault("molecular_weight", 0.0)
        nc.setdefault("tg", 0.0)
        nc.setdefault("pid", "")
        nc.setdefault("polymer_class", "")
        res["new_compound"] = nc

        # similar_compounds defaults
        scs = res.get("similar_compounds", [])
        for sc in scs:
            sc.setdefault("name", "")
            sc.setdefault("formula", "")
            sc.setdefault("molecular_weight", 0.0)
            
            # --- PERBAIKAN LOGIKA DISINI ---
            # Hapus baris ini: nc.setdefault("tg_justification", "") 
            # (Karena 'nc' tidak seharusnya diedit di dalam loop 'sc')
            
            sc.setdefault("tg", 0.0)
            sc.setdefault("pid", "")
            sc.setdefault("polymer_class", "")
            
            # Hitung persentase
            score = sc.get("similarity_score", 0.0)
            sc.setdefault("similarity_percent", float(score) * 100)
            
        res["similar_compounds"] = scs
        item["result"] = res

    return JSONResponse(content=jsonable_encoder(data))