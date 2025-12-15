import os
import re
import json
import time
from typing import List, Optional, Dict, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI

load_dotenv()

# =============================================================================
# UTIL
# =============================================================================
def _clean(text: str, max_len: int = 300) -> str:
    text = re.sub(r"\s+", " ", (text or "")).strip()
    return text[:max_len]


def _fix_degree_symbol(text: str) -> str:
    return (text or "").replace("Â°", "°").strip()


def _extract_json_obj(text: str) -> Optional[dict]:
    """
    Ambil objek JSON pertama dari output LLM.
    Handle codefence + teks tambahan.
    """
    if not text:
        return None

    t = text.strip()
    t = re.sub(r"^```(?:json)?\s*", "", t, flags=re.IGNORECASE)
    t = re.sub(r"\s*```$", "", t)

    start = t.find("{")
    end = t.rfind("}")
    if start == -1 or end == -1 or end <= start:
        return None

    cand = t[start : end + 1].strip()
    try:
        return json.loads(cand)
    except Exception:
        return None


def get_llm(
    model_name: str,
    timeout: int = 25,
    max_retries: int = 2,
    temperature: float = 0.3,
):
    api_key = os.getenv("GOOGLE_API_KEY")
    if not api_key:
        # ✅ fail fast biar gampang debug di deploy
        raise RuntimeError("GOOGLE_API_KEY belum diset di environment (Koyeb).")

    return ChatGoogleGenerativeAI(
        model=model_name,
        temperature=temperature,
        google_api_key=api_key,
        max_retries=max_retries,
        timeout=timeout,
    )


# =============================================================================
# LLM SINGLETONS
# =============================================================================
_llm_fast: Optional[ChatGoogleGenerativeAI] = None
_llm_tg: Optional[ChatGoogleGenerativeAI] = None


def llm_fast():
    """
    Untuk nama + justifikasi umum (boleh kreatif).
    """
    global _llm_fast
    if _llm_fast is None:
        _llm_fast = get_llm(
            model_name="gemini-2.5-flash",  # ✅ FIX: bukan 1.5-flash
            timeout=25,
            max_retries=2,
            temperature=0.7,
        )
    return _llm_fast


def llm_tg():
    """
    Khusus Tg: lebih deterministik.
    """
    global _llm_tg
    if _llm_tg is None:
        _llm_tg = get_llm(
            model_name="gemini-2.5-flash",  # ✅ FIX: bukan 1.5-flash
            timeout=30,
            max_retries=2,
            temperature=0.1,
        )
    return _llm_tg


def safe_invoke(prompt: str, fallback: str, max_len: int):
    try:
        resp = llm_fast().invoke(prompt)
        return _clean(getattr(resp, "content", ""), max_len=max_len)
    except Exception as e:
        print("LLM error:", repr(e))
        return fallback


def safe_invoke_tg(prompt: str, fallback: str, max_len: int):
    """
    Retry ringan biar gak kelamaan.
    """
    for attempt in range(2):
        try:
            resp = llm_tg().invoke(prompt)
            return _clean(getattr(resp, "content", ""), max_len=max_len)
        except Exception as e:
            print(f"LLM Tg attempt {attempt+1} error:", repr(e))
            time.sleep(1.2)
    return fallback


# =============================================================================
# SIMPLE IN-MEMORY CACHE (nama + justifikasi)
# =============================================================================
_llm_cache: Dict[str, Tuple[str, str]] = {}


def _get_cached(smiles: str) -> Optional[Tuple[str, str]]:
    return _llm_cache.get(smiles)


def _set_cached(smiles: str, name: str, justification: str):
    _llm_cache[smiles] = (name, justification)


# =============================================================================
# NAME + JUSTIFICATION
# =============================================================================
def generate_compound_name(smiles: str) -> str:
    smiles = (smiles or "").strip()
    cached = _get_cached(smiles)
    if cached and cached[0]:
        return cached[0]

    prompt = f"""Buat nama senyawa (IUPAC-like jika bisa) untuk SMILES ini: {smiles}
Aturan: maksimal 50 karakter, profesional, tanpa penjelasan.
Output: HANYA nama."""
    name = safe_invoke(prompt, fallback="GeneratedCompound", max_len=50)

    _set_cached(smiles, name, cached[1] if cached else "")
    return name


def generate_new_justification(smiles: str, compound_name: str) -> str:
    smiles = (smiles or "").strip()
    cached = _get_cached(smiles)
    if cached and cached[1]:
        return cached[1]

    prompt = f"""SMILES: {smiles}
Nama: {compound_name}
Tugas: Jelaskan 2 kalimat (Bahasa Indonesia) kenapa struktur ini unik/novel.
Output: HANYA justifikasi."""
    justif = safe_invoke(prompt, fallback="Analisis struktur sedang diproses.", max_len=600)

    _set_cached(smiles, compound_name, justif)
    return justif


def justify_similar_compounds_batch(compound_name: str, dataset_smiles_list: List[str]) -> List[str]:
    if not dataset_smiles_list:
        return []

    smiles_block = "\n".join([f"[{i+1}] {s}" for i, s in enumerate(dataset_smiles_list)])

    prompt = f"""Target: {compound_name}
Daftar SMILES:
{smiles_block}
Tugas: untuk tiap item, buat 1 kalimat Bahasa Indonesia kenapa mirip.
Format WAJIB:
[1] ...
[2] ...
[3] ..."""

    text = safe_invoke(prompt, fallback="", max_len=2000)

    results: List[str] = []
    for i in range(len(dataset_smiles_list)):
        m = re.search(rf"\[{i+1}\]\s*(.*?)(?=\[\d+\]|$)", text, flags=re.DOTALL)
        results.append(_clean(m.group(1), 300) if m else "Kemiripan struktur umum.")
    return results


# =============================================================================
# TG PREDICTION
# =============================================================================
def tg_fallback_heuristic(smiles: str) -> dict:
    """
    Fallback darurat (tanpa LLM) supaya tidak default 100 terus.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"tg": 0.0, "tg_justification": "Invalid SMILES"}

    mw = rdMolDescriptors.CalcExactMolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    rot = Lipinski.NumRotatableBonds(mol)

    tg = 30.0
    tg += 25.0 * rings
    tg -= 4.5 * rot
    if mw > 200:
        tg += 15.0

    tg = max(-100.0, min(300.0, tg))
    return {"tg": round(tg, 1), "tg_justification": "Fallback RDKit: estimasi ring/rotatable/MW."}


def predict_tg_with_llm(smiles: str, compound_name: str) -> dict:
    smiles = (smiles or "").strip()

    # 1) Diskcache hit
    try:
        from app.cache import cache, key_tg_pred
        cached = cache.get(key_tg_pred(smiles))
        if isinstance(cached, dict) and "tg" in cached:
            return {
                "tg": float(cached.get("tg", 0.0)),
                "tg_justification": _clean(str(cached.get("tg_justification", "")), 250),
            }
    except Exception:
        pass

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"tg": 0.0, "tg_justification": "SMILES error"}

    mw = rdMolDescriptors.CalcExactMolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)

    # Prompt ringkas = lebih cepat = lebih kecil timeout
    prompt = f"""Prediksi Tg (°C) untuk kandidat polimer berikut.
SMILES: {smiles}
Nama: {compound_name}
MW: {mw:.2f}
Rings: {rings}

Balas HANYA JSON valid:
{{"predicted_tg": -20.0, "tg_justification": "1-2 kalimat singkat."}}

Aturan:
- predicted_tg float -100..300
- tg_justification max 220 karakter
"""

    raw = safe_invoke_tg(prompt, fallback="", max_len=900)
    data = _extract_json_obj(raw)

    if isinstance(data, dict) and "predicted_tg" in data:
        try:
            tg = float(data.get("predicted_tg"))
            just = str(data.get("tg_justification", "AI Prediction"))
            tg = max(-100.0, min(300.0, tg))

            result = {
                "tg": round(tg, 1),
                "tg_justification": _clean(_fix_degree_symbol(just), 250),
            }

            try:
                from app.cache import cache, key_tg_pred
                cache.set(key_tg_pred(smiles), result, expire=60 * 60 * 24 * 7)
            except Exception:
                pass

            return result
        except Exception:
            pass

    # fallback
    return tg_fallback_heuristic(smiles)