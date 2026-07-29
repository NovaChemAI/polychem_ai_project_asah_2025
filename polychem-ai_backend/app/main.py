from contextlib import asynccontextmanager
import os
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Optional, Union

from app.schemas import RecommendRequest, RecommendResponse
from app.store import load_dataset
from app.core import recommend_new_compound
from app.settings import STATIC_DIR, COMPOUNDS_DIR
from app.cache import get_history
from app.auth_dependency import get_current_user, get_current_user_no_verify
from app.services.library import save_to_library
from app.services.library import remove_from_library
from app.auth_dependency import get_optional_user
from app.services.history import add_to_history
from app.services.library import get_user_library
from app.services.history import get_user_history
from app.services.users import sync_user_profile
from app.services.library import check_is_saved_by_smiles

import httpx
from app.settings import RECAPTCHA_SECRET_KEY

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

os.makedirs(STATIC_DIR, exist_ok=True)
os.makedirs(COMPOUNDS_DIR, exist_ok=True)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

frontend_url = os.getenv("FRONTEND_URL")
allow_origins = []
if frontend_url:
    allow_origins.append(frontend_url)

# dev origins
allow_origins += [
    "http://localhost:5173",
    "http://127.0.0.1:5173",
    "http://localhost:3000",
    "http://127.0.0.1:3000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# =========================================================
# Request Schemas
# =========================================================
class SaveLibraryRequest(BaseModel):
    id: Optional[str] = None
    name: Optional[str] = None
    smiles: Optional[str] = None
    category: Optional[str] = None
    image: Optional[str] = None
    image_url: Optional[str] = None
    properties: Optional[dict] = None
    score: Optional[Union[str, float]] = None
    isAiResult: Optional[bool] = False

class CaptchaRequest(BaseModel):
    captcha_token: str

# =========================================================
# API Endpoints
# =========================================================

@app.get("/health")
def health():
    return {"status": "ok"}

# =========================================================

@app.post("/library")
def add_to_library(req: SaveLibraryRequest, user=Depends(get_current_user)):
    uid = user["uid"]
    try:
        save_to_library(uid, req.dict())
        return {"success": True}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Gagal menyimpan: {str(e)}")

# =========================================================

@app.get("/library")
def list_library(user=Depends(get_current_user)):
    uid = user["uid"]
    try:
        items = get_user_library(uid)
        return JSONResponse(content=jsonable_encoder(items))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Gagal mengambil data: {str(e)}")

# =========================================================

@app.delete("/library/{item_id}")
def delete_from_library(item_id: str, user=Depends(get_current_user)):
    uid = user["uid"]
    try:
        remove_from_library(uid, item_id)
        return {"success": True}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Gagal menghapus: {str(e)}")

# =========================================================

@app.get("/library/check")
def check_saved_status(smiles: str, user=Depends(get_current_user)):
    uid = user["uid"]
    try:
        result = check_is_saved_by_smiles(uid, smiles)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Gagal cek status: {str(e)}")
    
# =========================================================
# NOTE: Hanya ADA SATU definisi /predict di sini.
# Sebelumnya ada dua @app.post("/predict") yang menyebabkan
# definisi pertama (tanpa penyimpanan history) selalu dipakai,
# sehingga add_to_history() tidak pernah terpanggil.
# =========================================================

@app.post("/predict", response_model=RecommendResponse)
def predict(req: RecommendRequest, user=Depends(get_optional_user)):
    try:
        result = recommend_new_compound(req.smiles)
    except Exception as e:
        print(f"Error di pipeline: {e}")
        raise HTTPException(status_code=503, detail=f"Pipeline error: {str(e)}")

    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])

    if user:
        try:
            add_to_history(user["uid"], req.smiles, result)
        except Exception as e:
            print(f"⚠️ Gagal simpan history: {e}")

    return result

# =========================================================

@app.get("/history/mine")
def my_history(user=Depends(get_current_user)):
    uid = user["uid"]
    try:
        items = get_user_history(uid)
        return JSONResponse(content=jsonable_encoder(items))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Gagal mengambil riwayat: {str(e)}")
    
# =========================================================

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
    
# =========================================================

class SyncProfileRequest(BaseModel):
    name: Optional[str] = None
    photo_url: Optional[str] = None

@app.post("/auth/sync-profile")
def sync_profile(req: SyncProfileRequest, user=Depends(get_current_user_no_verify)):
    uid = user["uid"]
    email = user.get("email", "")
    name = req.name or user.get("name", "Unknown")
    try:
        profile = sync_user_profile(uid, name, email, req.photo_url or "")
        return {"success": True, "profile": profile}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Gagal sinkronisasi profil: {str(e)}")
    
# =========================================================

@app.post("/auth/verify-captcha")
async def verify_captcha(req: CaptchaRequest):
    if not RECAPTCHA_SECRET_KEY:
        raise HTTPException(status_code=500, detail="RECAPTCHA_SECRET_KEY belum dikonfigurasi di server")

    async with httpx.AsyncClient() as client:
        resp = await client.post(
            "https://www.google.com/recaptcha/api/siteverify",
            data={"secret": RECAPTCHA_SECRET_KEY, "response": req.captcha_token},
        )
        result = resp.json()

    if not result.get("success"):
        raise HTTPException(status_code=400, detail="Verifikasi CAPTCHA gagal")

    return {"success": True}