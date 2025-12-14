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


def get_llm(
    model_name: str,
    timeout: int = 60,
    max_retries: int = 3,
    temperature: float = 0.3
):
    api_key = os.getenv("GOOGLE_API_KEY")
    if not api_key:
        print("⚠️ WARNING: GOOGLE_API_KEY belum diset.")

    return ChatGoogleGenerativeAI(
        model=model_name,
        temperature=temperature,
        google_api_key=api_key,
        max_retries=max_retries,
        timeout=timeout,
        request_timeout=timeout
    )


# =============================================================================
# LLM SINGLETONS
# =============================================================================
_llm_fast: Optional[ChatGoogleGenerativeAI] = None
_llm_tg: Optional[ChatGoogleGenerativeAI] = None


def llm_fast():
    global _llm_fast
    if _llm_fast is None:
        _llm_fast = get_llm(
            "gemini-1.5-flash", 
            timeout=40,
            temperature=0.7
        )
    return _llm_fast


def llm_tg():
    global _llm_tg
    if _llm_tg is None:
        _llm_tg = get_llm(
            "gemini-1.5-flash", 
            timeout=60, 
            temperature=0.1 
        )
    return _llm_tg


def safe_invoke(prompt: str, fallback: str, max_len: int):
    try:
        resp = llm_fast().invoke(prompt)
        return _clean(getattr(resp, "content", ""), max_len=max_len)
    except Exception as e:
        print(f"❌ LLM Name Error: {e}")
        return fallback


def safe_invoke_tg(prompt: str, fallback: str, max_len: int):
    # Retry manual logic
    for attempt in range(2):
        try:
            resp = llm_tg().invoke(prompt)
            return _clean(getattr(resp, "content", ""), max_len=max_len)
        except Exception as e:
            print(f"⚠️ Percobaan Tg ke-{attempt+1} gagal: {e}")
            time.sleep(2) 
            
    print("❌ LLM Tg Gagal Total.")
    return fallback


# =============================================================================
# CACHE
# =============================================================================
_llm_cache: Dict[str, Tuple[str, str]] = {}

def _get_cached(smiles: str) -> Optional[Tuple[str, str]]:
    return _llm_cache.get(smiles)

def _set_cached(smiles: str, name: str, justification: str):
    _llm_cache[smiles] = (name, justification)


# =============================================================================
# GENERATORS
# =============================================================================
def generate_compound_name(smiles: str) -> str:
    smiles = (smiles or "").strip()
    cached = _get_cached(smiles)
    if cached and cached[0]: return cached[0]

    prompt = f"""Task: Generate a scientific IUPAC-like name for this SMILES: {smiles}
Constraint: Max 50 chars, Unique, Professional.
Output: ONLY the name."""
    
    name = safe_invoke(prompt, fallback="GeneratedCompound", max_len=50)
    _set_cached(smiles, name, cached[1] if cached else "")
    return name


def generate_new_justification(smiles: str, compound_name: str) -> str:
    smiles = (smiles or "").strip()
    cached = _get_cached(smiles)
    if cached and cached[1]: return cached[1]

    prompt = f"""Structure: {smiles} ({compound_name})
Task: Explain in 2 sentences (Indonesian) why this compound is unique/novel.
Focus on: Functional groups and potential stability.
Output: Justification text only."""

    justif = safe_invoke(prompt, fallback="Analisis struktur sedang diproses.", max_len=600)
    _set_cached(smiles, compound_name, justif)
    return justif


def justify_similar_compounds_batch(compound_name: str, dataset_smiles_list: List[str]) -> List[str]:
    if not dataset_smiles_list: return []
    smiles_block = "\n".join([f"[{i+1}] {s}" for i, s in enumerate(dataset_smiles_list)])
    prompt = f"""Target: {compound_name}
List:
{smiles_block}
Task: Give 1 sentence (Indonesian) for each, explaining structural similarity.
Format: [1] text [2] text"""
    
    text = safe_invoke(prompt, fallback="", max_len=2000)
    results = []
    for i in range(len(dataset_smiles_list)):
        match = re.search(rf"\[{i+1}\]\s*(.*?)(?=\[\d+\]|$)", text, re.DOTALL)
        results.append(_clean(match.group(1), 300) if match else "Kemiripan struktur umum.")
    return results


# =============================================================================
# TG PREDICTION
# =============================================================================

def _extract_json_obj(text: str) -> Optional[dict]:
    if not text: return None
    t = text.strip()
    t = re.sub(r"^```(?:json)?\s*", "", t, flags=re.IGNORECASE)
    t = re.sub(r"\s*```$", "", t)
    start = t.find("{")
    end = t.rfind("}")
    if start == -1 or end == -1: return None
    try: return json.loads(t[start:end+1])
    except: return None


# --- PERBAIKAN UTAMA: HAPUS UNDERSCORE (_) AGAR BISA DIIMPORT ---
def tg_fallback_heuristic(smiles: str) -> dict:
    """
    Fallback darurat hanya jika AI mati total.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return {"tg": 0.0, "tg_justification": "Invalid SMILES"}
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    
    # Rumus Heuristik
    tg = 50.0 
    tg += 30.0 * rings 
    tg -= 5.0 * rot    
    if mw > 200: tg += 20
    
    return {
        "tg": round(tg, 1),
        "tg_justification": "Estimasi kalkulasi (AI Timeout)."
    }


def predict_tg_with_llm(smiles: str, compound_name: str) -> dict:
    """
    Fungsi Utama Prediksi Tg.
    """
    smiles = (smiles or "").strip()

    # 1. CEK CACHE
    try:
        from app.cache import cache, key_tg_pred
        cached = cache.get(key_tg_pred(smiles))
        if cached and cached.get("tg", 0) != 0:
            print(f"✅ Cache Hit: {cached['tg']}")
            return cached
    except: pass

    mol = Chem.MolFromSmiles(smiles)
    if not mol: return {"tg": 0.0, "tg_justification": "SMILES error"}

    mw = rdMolDescriptors.CalcExactMolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    
    print(f"🚀 AI Calculating Tg for: {smiles[:15]}...")

    # 2. PROMPT
    prompt = f"""Role: Expert Polymer Physicist.
Task: Predict the Glass Transition Temperature (Tg) in Celsius.

Compound Info:
- SMILES: {smiles}
- Name: {compound_name}
- MW: {mw:.2f}
- Rings: {rings}

Analysis Rules:
1. Consider chain stiffness (aromatic rings = high Tg).
2. Consider intermolecular forces.
3. If rigid polymer, Tg usually > 80°C.
4. If flexible rubber, Tg < 0°C.

Output ONLY JSON:
{{
  "predicted_tg": <float_value>,
  "tg_justification": "<Explain strictly based on structure rigidity/flexibility>"
}}
"""

    raw = safe_invoke_tg(prompt, fallback="", max_len=1000)
    data = _extract_json_obj(raw)

    if data and "predicted_tg" in data:
        try:
            tg = float(data["predicted_tg"])
            just = data.get("tg_justification", "AI Prediction")
            tg = max(-100.0, min(400.0, tg))
            
            result = {
                "tg": round(tg, 1),
                "tg_justification": _clean(_fix_degree_symbol(just), 250)
            }
            
            try:
                from app.cache import cache, key_tg_pred
                cache.set(key_tg_pred(smiles), result, expire=60*60*24*7)
            except: pass
            
            print(f"✅ AI Success: {result['tg']}°C")
            return result
        except: pass

    print("⚠️ AI Failed. Using Fallback.")
    # Panggil fungsi yang sekarang sudah PUBLIK
    return tg_fallback_heuristic(smiles)