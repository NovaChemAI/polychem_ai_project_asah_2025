import os
import re
from typing import List, Optional, Dict, Tuple

from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI

load_dotenv()

# =============================================================================
# UTIL
# =============================================================================
def _clean(text: str, max_len: int = 300) -> str:
    text = re.sub(r"\s+", " ", text).strip()
    return text[:max_len]

def get_llm(model_name: str, timeout: int = 12):
    api_key = os.getenv("GOOGLE_API_KEY")
    if not api_key:
        raise RuntimeError(
            "GOOGLE_API_KEY belum diset. Isi di file .env atau environment."
        )

    return ChatGoogleGenerativeAI(
        model=model_name,
        temperature=0.7,
        google_api_key=api_key,
        max_retries=0, 
        timeout=timeout, 
    )

# LLM SINGLETON (Tier 1)
_llm_fast: Optional[ChatGoogleGenerativeAI] = None

def llm_fast():
    global _llm_fast
    if _llm_fast is None:
        _llm_fast = get_llm("gemini-2.5-flash", timeout=12)
    return _llm_fast

def safe_invoke(prompt: str, fallback: str, max_len: int):
    try:
        resp = llm_fast().invoke(prompt)
        return _clean(resp.content, max_len=max_len)
    except Exception:
        return fallback

# SIMPLE IN-MEMORY CACHE
_llm_cache: Dict[str, Tuple[str, str]] = {}

def _get_cached(smiles: str) -> Optional[Tuple[str, str]]:
    return _llm_cache.get(smiles)

def _set_cached(smiles: str, name: str, justification: str):
    _llm_cache[smiles] = (name, justification)

# PUBLIC API
def generate_compound_name(smiles: str) -> str:
    smiles = smiles.strip()

    # cache hit
    cached = _get_cached(smiles)
    if cached:
        name, _ = cached
        if name:
            return name

    prompt = f"""Buat nama senyawa kimia baru yang belum teridentifikasi untuk SMILES ini: {smiles}

Nama harus:
1. Unik dan belum ada di dunia
2. Mengikuti konvensi IUPAC (jika memungkinkan)
3. Menarik dan mudah diingat
4. Maksimal 50 karakter

Hanya berikan nama, tanpa penjelasan tambahan."""

    name = safe_invoke(prompt, fallback="DummyName", max_len=50)

    _set_cached(smiles, name, cached[1] if cached else "")
    return name


def generate_new_justification(smiles: str, compound_name: str) -> str:
    smiles = smiles.strip()

    # cache hit
    cached = _get_cached(smiles)
    if cached:
        _, justif = cached
        if justif:
            return justif

    prompt = f"""Berdasarkan struktur SMILES: {smiles}
Nama senyawa: {compound_name}

Berikan justifikasi singkat (2-3 kalimat dalam Bahasa Indonesia) mengapa senyawa ini unik dan belum teridentifikasi.
Pertimbangkan: struktur molekul, potensi aplikasi, dan keunikan karakteristik kimia.

Hanya berikan justifikasi tanpa nomor atau penjelasan tambahan."""

    justif = safe_invoke(
        prompt,
        fallback="Justifikasi sementara (LLM timeout/off).",
        max_len=500
    )

    _set_cached(smiles, compound_name, justif)
    return justif


def justify_similar_compounds_batch(
    compound_name: str,
    dataset_smiles_list: List[str],
) -> List[str]:
    if not dataset_smiles_list:
        return []

    smiles_block = "\n".join(
        [f"[{i+1}] {s}" for i, s in enumerate(dataset_smiles_list)]
    )

    prompt = f"""Target compound: {compound_name}

Dataset compound SMILES:
{smiles_block}

Tugas:
Untuk setiap SMILES di atas, berikan justifikasi singkat (2-3 kalimat Bahasa Indonesia) mengapa mirip dengan target.

Aturan output:
- Tulis persis dalam format:
[1] <justifikasi>
[2] <justifikasi>
[3] <justifikasi>
- Jangan tambah nomor lain, jangan tambah teks di luar format itu."""

    text = safe_invoke(prompt, fallback="", max_len=2000)

    results: List[str] = []
    for i in range(len(dataset_smiles_list)):
        pattern = rf"\[{i+1}\]\s*(.*?)(?=\[\d+\]|$)"
        match = re.search(pattern, text, flags=re.DOTALL)
        if match:
            results.append(_clean(match.group(1), max_len=400))
        else:
            results.append("Mirip secara struktur umum (fallback).")

    return results
