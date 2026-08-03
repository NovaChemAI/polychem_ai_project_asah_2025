import os
import re
import json
import time
from typing import List, Optional, Dict, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski, Descriptors
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI

load_dotenv()

# =============================================================================
# UTIL
# =============================================================================
def _clean(text: str, max_len: int = 300) -> str:
    """
    Bersihkan whitespace, lalu potong ke max_len TANPA memotong di tengah kata
    atau kalimat kalau bisa dihindari.
    """
    text = re.sub(r"\s+", " ", (text or "")).strip()
    if len(text) <= max_len:
        return text

    truncated = text[:max_len]

    # Coba potong mundur ke akhir kalimat/klausa terdekat (". " atau ", ")
    last_period = max(truncated.rfind(". "), truncated.rfind(", "))
    if last_period > max_len * 0.6:
        return truncated[: last_period + 1].strip()

    # Kalau tidak ada, potong mundur ke spasi terakhir supaya tidak motong kata
    last_space = truncated.rfind(" ")
    if last_space > max_len * 0.6:
        return truncated[:last_space].strip() + "..."

    return truncated.strip() + "..."


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
        # coba ambil blok JSON paling "lebar" kalau ada noise
        m = re.search(r"\{.*\}", t, flags=re.DOTALL)
        if not m:
            return None
        try:
            return json.loads(m.group(0))
        except Exception:
            return None


def get_llm(
    model_name: str,
    timeout: int = 45,
    max_retries: int = 2,
    temperature: float = 0.0,
):
    api_key = os.getenv("GOOGLE_API_KEY") or os.getenv("GEMINI_API_KEY")

    if not api_key:
        raise RuntimeError("GOOGLE_API_KEY/GEMINI_API_KEY belum diset di environment.")

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
    Untuk nama + justifikasi umum.
    Model bisa dioverride via ENV: GEMINI_MODEL_FAST
    """
    global _llm_fast
    if _llm_fast is None:
        model = os.getenv("GEMINI_MODEL_FAST", "gemini-2.5-flash")
        _llm_fast = get_llm(
            model_name=model,
            timeout=45,
            max_retries=2,
            temperature=0.0,
        )
    return _llm_fast


def llm_tg():
    """
    Khusus Tg: deterministik.
    Model bisa dioverride via ENV: GEMINI_MODEL_TG
    """
    global _llm_tg
    if _llm_tg is None:
        model = os.getenv("GEMINI_MODEL_TG", "gemini-2.5-flash")
        _llm_tg = get_llm(
            model_name=model,
            timeout=45,
            max_retries=2,
            temperature=0.0,
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
    try:
        resp = llm_tg().invoke(prompt)
        return _clean(getattr(resp, "content", ""), max_len=max_len)
    except Exception as e:
        print("LLM Tg error:", repr(e))
        return fallback


# =============================================================================
# CACHE (nama + justifikasi) — pakai diskcache PERSISTEN
# =============================================================================
NAME_JUSTIF_TTL_SECONDS = 60 * 60 * 24 * 7  # 7 hari


def _name_justif_key(smiles: str) -> str:
    try:
        from app.cache import CACHE_VERSION
    except Exception:
        CACHE_VERSION = "v2"
    smiles_n = " ".join((smiles or "").strip().split())
    return f"{CACHE_VERSION}::name_justif::{smiles_n}"


def _get_cached(smiles: str) -> Optional[Tuple[str, str]]:
    try:
        from app.cache import cache
        data = cache.get(_name_justif_key(smiles))
        if isinstance(data, (list, tuple)) and len(data) == 2 and data[0] and data[1]:
            return (str(data[0]), str(data[1]))
    except Exception as e:
        print("⚠️ _get_cached error:", repr(e))
    return None


def _set_cached(smiles: str, name: str, justification: str):
    try:
        from app.cache import cache
        cache.set(
            _name_justif_key(smiles),
            (name, justification),
            expire=NAME_JUSTIF_TTL_SECONDS,
        )
    except Exception as e:
        print("⚠️ _set_cached error:", repr(e))


# =============================================================================
# FALLBACK JUSTIFICATION PANJANG (RDKit)
# =============================================================================
def fallback_justification_long(smiles: str, compound_name: str = "") -> str:
    """
    Fallback panjang berbasis RDKit supaya output tetap informatif walau LLM gagal/timeout.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return (
            "Struktur SMILES tidak valid sehingga analisis otomatis tidak dapat dilakukan. "
            "Periksa kembali format SMILES (terutama karakter khusus seperti '*' atau penulisan ikatan '='). "
            "Jika SMILES sudah valid, sistem akan menurunkan fitur struktur dan memberikan justifikasi yang lebih detail."
        )

    mw = Descriptors.MolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    name_part = f" ({compound_name})" if compound_name else ""
    return (
        f"Struktur{name_part} dapat dijelaskan melalui keseimbangan antara kekakuan segmen (jumlah cincin {rings}) "
        f"dan fleksibilitas rantai (ikatan rotatable {rot}), yang berpengaruh langsung terhadap mobilitas rantai pada material. "
        f"Dengan massa molekul sekitar {mw:.1f} g/mol dan polaritas yang tercermin dari HBD={hbd} serta HBA={hba}, "
        f"senyawa ini berpotensi menunjukkan interaksi antarmolekul yang memengaruhi sifat seperti Tg, kompatibilitas, dan kelarutan. "
        f"Kombinasi motif struktur, tingkat rotasi ikatan, dan gugus fungsi menjadikannya kandidat yang bisa berbeda dari analog sederhana "
        f"pada perilaku fisik maupun potensi aplikasinya."
    )


# =============================================================================
# NAME + JUSTIFICATION (DIGABUNG JADI 1 LLM CALL)
# =============================================================================
def generate_name_and_justification(smiles: str) -> Tuple[str, str]:
    """
    Satu panggilan LLM untuk nama + justifikasi sekaligus (efisiensi).
    Hasil di-cache secara persisten (diskcache) per SMILES.
    """
    smiles = (smiles or "").strip()

    cached = _get_cached(smiles)
    if cached:
        return cached

    prompt = f"""SMILES: {smiles}

Tugas:
1. Buat nama senyawa (IUPAC-like jika bisa), maksimal 50 karakter, profesional, tanpa penjelasan.
2. Tulis justifikasi 3-6 kalimat (Bahasa Indonesia) mengapa struktur ini unik/novel.
   WAJIB menyebut minimal 2 aspek berikut:
   - gugus fungsi / polaritas
   - fleksibilitas rantai (rotatable bonds)
   - kekakuan (cincin/aromatik)
   - implikasi terhadap sifat material/polimer (misal Tg, kompatibilitas, stabilitas)

Balas HANYA JSON valid, tanpa markdown/codefence, dengan format persis:
{{"name": "...", "justification": "..."}}
"""

    fallback_justif = fallback_justification_long(smiles, "")
    raw = safe_invoke(prompt, fallback="", max_len=1600)
    data = _extract_json_obj(raw)

    if isinstance(data, dict) and data.get("name") and data.get("justification"):
        name = _clean(str(data["name"]), 50)
        justif = _clean(str(data["justification"]), 1200)
        if len(justif) < 220:
            justif = _clean(fallback_justif, 1200)
    else:
        name = "GeneratedCompound"
        justif = fallback_justif

    _set_cached(smiles, name, justif)
    return name, justif


def generate_compound_name(smiles: str) -> str:
    name, _ = generate_name_and_justification(smiles)
    return name


def generate_new_justification(smiles: str, compound_name: str) -> str:
    _, justif = generate_name_and_justification(smiles)
    return justif


# =============================================================================
# JUSTIFIKASI SENYAWA MIRIP (BATCH)
# =============================================================================
def justify_similar_compounds_batch(compound_name: str, dataset_smiles_list: List[str]) -> List[str]:
    if not dataset_smiles_list:
        return []

    smiles_block = "\n".join([f"[{i+1}] {s}" for i, s in enumerate(dataset_smiles_list)])

    prompt = f"""Target: {compound_name}

Daftar SMILES:
{smiles_block}

Tugas:
Untuk setiap SMILES, tulis justifikasi 2-3 kalimat (Bahasa Indonesia) mengapa mirip dengan target.
Sebutkan motif yang sama bila ada: ester/eter/aromatik/percabangan/polaritas, atau kemiripan pola rantai.

WAJIB ikuti format ini PERSIS, satu baris per nomor, tanpa teks pembuka/penutup lain:
[1] <justifikasi>
[2] <justifikasi>
[3] <justifikasi>"""

    text = safe_invoke(prompt, fallback="", max_len=5000)

    # --- DEBUG: cek jawaban mentah kalau masih ada masalah parsing ---
    # print("=== RAW GEMINI (similar compounds) ===")
    # print(repr(text))

    results: List[str] = _parse_numbered_justifs(text, dataset_smiles_list, compound_name)
    return results


def _parse_numbered_justifs(text: str, dataset_smiles_list: List[str], compound_name: str) -> List[str]:
    """
    Parsing yang lebih toleran:
    1. Coba format [1]/[2]/[3] (regex asli).
    2. Kalau gagal, coba format "1." / "1)" (variasi umum lain dari Gemini).
    3. Kalau masih gagal tapi teks cukup panjang & tidak ada nomor sama sekali,
       coba split per baris sebagai upaya terakhir sebelum fallback.
    4. Ambang minimal panjang diturunkan (200 -> 60) supaya jawaban valid
       yang ringkas tidak ikut dibuang.
    """
    n = len(dataset_smiles_list)
    results: List[Optional[str]] = []

    if not text or not text.strip():
        return [_generic_fallback(compound_name) for _ in range(n)]

    patterns = [
        r"\[{i}\]\s*(.*?)(?=\[\d+\]|$)",                              # [1] ...
        r"(?:^|\n)\s*{i}[\.\)]\s*(.*?)(?=(?:\n\s*\d+[\.\)])|$)",       # 1. ... atau 1) ...
    ]

    for i in range(n):
        cand = None
        for pat in patterns:
            m = re.search(pat.format(i=i + 1), text, flags=re.DOTALL)
            if m:
                candidate_text = _clean(m.group(1), 900)
                if len(candidate_text) >= 60:
                    cand = candidate_text
                    break
        results.append(cand)

    if all(r is None for r in results):
        lines = [l.strip() for l in text.split("\n") if len(l.strip()) >= 60]
        for i in range(n):
            if i < len(lines):
                results[i] = _clean(lines[i], 900)

    final_results: List[str] = []
    for i in range(n):
        final_results.append(results[i] if results[i] else _generic_fallback(compound_name))

    return final_results


def _generic_fallback(compound_name: str) -> str:
    return _clean(
        f"Senyawa ini menunjukkan kemiripan dengan {compound_name} melalui pola ikatan dan motif gugus fungsi yang sejenis. "
        f"Kedekatan fitur seperti panjang rantai, percabangan, atau segmen polar dapat menghasilkan fingerprint yang berdekatan. "
        f"Kesamaan tersebut sering berkorelasi dengan kecenderungan sifat fisik yang mirip pada level segmen rantai.",
        900,
    )


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

    prompt = f"""Prediksi Tg (°C) untuk kandidat polimer berikut.
SMILES: {smiles}
Nama: {compound_name}
MW: {mw:.2f}
Rings: {rings}

Balas HANYA JSON valid (tg_justification WAJIB dalam Bahasa Indonesia):
{{"predicted_tg": -20.0, "tg_justification": "2-4 kalimat alasan berbasis kekakuan/fleksibilitas rantai dan polaritas, dalam Bahasa Indonesia."}}

Aturan:
- predicted_tg float -100..300
- tg_justification maksimal 480 karakter, Bahasa Indonesia
"""

    raw = safe_invoke_tg(prompt, fallback="", max_len=1400)
    data = _extract_json_obj(raw)

    if isinstance(data, dict) and "predicted_tg" in data:
        try:
            tg = float(data.get("predicted_tg"))
            just = str(data.get("tg_justification", "AI Prediction"))
            tg = max(-100.0, min(300.0, tg))

            result = {
                "tg": round(tg, 1),
                "tg_justification": _clean(_fix_degree_symbol(just), 480),
            }

            try:
                from app.cache import cache, key_tg_pred
                cache.set(key_tg_pred(smiles), result, expire=60 * 60 * 24 * 7)
            except Exception:
                pass

            return result
        except Exception:
            pass

    return tg_fallback_heuristic(smiles)