import os
import time
import sqlite3
from diskcache import Cache
from typing import List, Dict, Any, Optional

# =========================
# CACHE DIR (writable in Koyeb)
# =========================
def _ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path

def _pick_cache_dir() -> str:
    """
    Default pakai /tmp karena di Koyeb folder project (/app) bisa read-only.
    Bisa override lewat env CACHE_DIR.
    """
    env_dir = os.getenv("CACHE_DIR")
    if env_dir:
        return _ensure_dir(env_dir)

    return _ensure_dir("/tmp/polychem_cache")

CACHE_DIR = _pick_cache_dir()

# bump ini kalau prompt/model/dataset berubah besar
CACHE_VERSION = os.getenv("CACHE_VERSION", "v1")

# history config - optimized for Nano instance
HISTORY_KEY = f"{CACHE_VERSION}::history"
HISTORY_LIMIT = int(os.getenv("HISTORY_LIMIT", "5"))  # Reduced from 10 to 5 for memory
HISTORY_TTL_SECONDS = int(os.getenv("HISTORY_TTL_SECONDS", str(60 * 60)))  # 1 jam

# =========================
# Cache object global
# =========================
def _open_cache(cache_dir: str) -> Cache:
    """
    Buka diskcache. Kalau kena readonly sqlite, fallback ke /tmp.
    
    SIZE_LIMIT: 50MB (reduced from 300MB for Nano instance memory optimization)
    """
    try:
        cache_obj = Cache(cache_dir, size_limit=int(5e7))  # 50MB
        # Cleanup old entries on startup to free memory
        try:
            cache_obj.evict(ratio=0.3)  # Remove 30% oldest entries if at limit
            print(f"✅ Cache cleanup done. Cache dir: {cache_dir}")
        except Exception:
            pass
        return cache_obj
    except sqlite3.OperationalError as e:
        # fallback keras kalau ternyata readonly
        fallback_dir = _ensure_dir("/tmp/polychem_cache_fallback")
        print(f"⚠️ diskcache readonly at {cache_dir}. Fallback to {fallback_dir}. Error: {repr(e)}")
        cache_obj = Cache(fallback_dir, size_limit=int(5e7))  # 50MB
        try:
            cache_obj.evict(ratio=0.3)
        except Exception:
            pass
        return cache_obj

cache = _open_cache(CACHE_DIR)

# =========================
# HELPERS
# =========================
def _norm(s: Optional[str]) -> str:
    """Normalize biar key stabil (hapus spasi berlebih)."""
    return " ".join(str(s or "").strip().split())

def key_new_compound(smiles: str) -> str:
    smiles_n = _norm(smiles)
    return f"{CACHE_VERSION}::new::{smiles_n}"

def key_similar_justif(compound_name: str, top_smiles: List[str]) -> str:
    cname_n = _norm(compound_name)
    joined = "|".join(_norm(s) for s in (top_smiles or []))
    return f"{CACHE_VERSION}::similar::{cname_n}::{joined}"

def key_tg_pred(smiles: str) -> str:
    smiles_n = _norm(smiles)
    return f"{CACHE_VERSION}::tg::{smiles_n}"

# =========================
# HISTORY (TTL-based)
# =========================
def _prune_history(hist: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    now = time.time()
    fresh: List[Dict[str, Any]] = []
    for item in hist:
        ts = item.get("ts")
        if isinstance(ts, (int, float)) and (now - ts) <= HISTORY_TTL_SECONDS:
            fresh.append(item)
    return fresh

def push_history(item: Dict[str, Any]) -> None:
    """
    Simpan history request terbaru (aman: cache error tidak bikin API crash).
    """
    try:
        hist = cache.get(HISTORY_KEY, default=[])
        if not isinstance(hist, list):
            hist = []

        hist = _prune_history(hist)

        new_item = dict(item)
        new_item["ts"] = time.time()

        hist.insert(0, new_item)
        hist = hist[:HISTORY_LIMIT]

        cache.set(HISTORY_KEY, hist, expire=HISTORY_TTL_SECONDS)
    except Exception as e:
        print("⚠️ push_history cache error:", repr(e))

def get_history() -> List[Dict[str, Any]]:
    """
    Ambil list history (aman: cache error -> return []).
    """
    try:
        hist = cache.get(HISTORY_KEY, default=[])
        if not isinstance(hist, list):
            return []

        hist2 = _prune_history(hist)
        if len(hist2) != len(hist):
            cache.set(HISTORY_KEY, hist2, expire=HISTORY_TTL_SECONDS)

        return hist2
    except Exception as e:
        print("⚠️ get_history cache error:", repr(e))
        return []