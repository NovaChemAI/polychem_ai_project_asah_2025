import os
import time
from diskcache import Cache
from typing import List, Dict, Any

# root project
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

# folder cache di root project: ./cache_data
CACHE_DIR = os.path.join(BASE_DIR, "cache_data")
os.makedirs(CACHE_DIR, exist_ok=True)

# bump ini kalau prompt/model/dataset berubah besar
CACHE_VERSION = "v1"

# history config
HISTORY_KEY = f"{CACHE_VERSION}::history"
HISTORY_LIMIT = 10
HISTORY_TTL_SECONDS = 60 * 60 # 1 jam

# Cache object global (300mb)
cache = Cache(CACHE_DIR, size_limit=int(3e8))

# HELPERS
def _norm(s: str) -> str:
    """Normalize biar key stabil (hapus spasi berlebih)."""
    return " ".join(s.strip().split())


def key_new_compound(smiles: str) -> str:
    """Key untuk cache nama+justifikasi senyawa baru."""
    smiles_n = _norm(smiles)
    return f"{CACHE_VERSION}::new::{smiles_n}"


def key_similar_justif(compound_name: str, top_smiles: List[str]) -> str:
    """Key untuk cache batch justifikasi similar compounds."""
    cname_n = _norm(compound_name)
    joined = "|".join(_norm(s) for s in top_smiles)
    return f"{CACHE_VERSION}::similar::{cname_n}::{joined}"

# HISTORY (TTL-based)
def _prune_history(hist: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Buang item yang sudah lewat TTL."""
    now = time.time()
    fresh = []
    for item in hist:
        ts = item.get("ts")
        if isinstance(ts, (int, float)) and (now - ts) <= HISTORY_TTL_SECONDS:
            fresh.append(item)
    return fresh


def push_history(item: Dict[str, Any]) -> None:
    """
    Simpan history request terbaru:
    - item diberi timestamp ts
    - prune berdasarkan TTL
    - limit 50
    - key history juga diset expire (TTL) supaya auto hilang kalau tidak update sama sekali
    """
    hist = cache.get(HISTORY_KEY, default=[])
    if not isinstance(hist, list):
        hist = []

    hist = _prune_history(hist)
    item = dict(item)
    item["ts"] = time.time()

    hist.insert(0, item)    
    hist = hist[:HISTORY_LIMIT]

    #expire untuk seluruh history key
    cache.set(HISTORY_KEY, hist, expire=HISTORY_TTL_SECONDS)


def get_history() -> List[Dict[str, Any]]:
    """
    Ambil list history:
    - prune TTL saat dibaca juga
    - update cache biar item expired beneran kebuang
    """
    hist = cache.get(HISTORY_KEY, default=[])
    if not isinstance(hist, list):
        return []

    hist2 = _prune_history(hist)

    if len(hist2) != len(hist):
        cache.set(HISTORY_KEY, hist2, expire=HISTORY_TTL_SECONDS)

    return hist2

def key_tg_pred(smiles: str) -> str:
    smiles_n = _norm(smiles)
    return f"{CACHE_VERSION}::tg::{smiles_n}"