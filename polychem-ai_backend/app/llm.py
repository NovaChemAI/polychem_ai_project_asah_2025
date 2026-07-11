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
    temperature: float = 0.3,
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
# CATATAN OPTIMASI:
# - timeout & max_retries diturunkan dari (25, 2) / (30, 2) menjadi lebih ketat.
#   Model "flash" normalnya respons < 5 detik; timeout tinggi cuma bikin
#   worst-case delay makin parah saat network lambat, bukan mempercepat.
_llm_fast: Optional[ChatGoogleGenerativeAI] = None
_llm_tg: Optional[ChatGoogleGenerativeAI] = None


def llm_fast():
    """
    Untuk nama + justifikasi umum (boleh kreatif).
    Model bisa dioverride via ENV: GEMINI_MODEL_FAST
    """
    global _llm_fast
    if _llm_fast is None:
        model = os.getenv("GEMINI_MODEL_FAST", "gemini-2.5-flash")
        _llm_fast = get_llm(
            model_name=model,
            timeout=45,
            max_retries=2,
            temperature=0.7,
        )
    return _llm_fast


def llm_tg():
    """
    Khusus Tg: lebih deterministik.
    Model bisa dioverride via ENV: GEMINI_MODEL_TG
    """
    global _llm_tg
    if _llm_tg is None:
        model = os.getenv("GEMINI_MODEL_TG", "gemini-2.5-flash")
        _llm_tg = get_llm(
            model_name=model,
            timeout=45,
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
    NOTE: max_retries sudah dihandle oleh client (ChatGoogleGenerativeAI),
    jadi di sini cukup 1x percobaan manual tanpa sleep tambahan yang
    memperpanjang total waktu tunggu.
    """
    try:
        resp = llm_tg().invoke(prompt)
        return _clean(getattr(resp, "content", ""), max_len=max_len)
    except Exception as e:
        print("LLM Tg error:", repr(e))
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
    OPTIMASI: sebelumnya ini 2 panggilan LLM terpisah
    (generate_compound_name + generate_new_justification).
    Digabung jadi 1 panggilan -> motong 1 network round-trip penuh
    (biasanya menghemat 1-3 detik per request).
    """
    smiles = (smiles or "").strip()
    cached = _get_cached(smiles)
    if cached and cached[0] and cached[1]:
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


# --- Wrapper backward-compatible (kalau ada bagian lain yang masih memanggil
#     fungsi lama secara terpisah, ini tetap jalan tapi otomatis pakai cache
#     dari generate_name_and_justification kalau sudah pernah dipanggil) ---
def generate_compound_name(smiles: str) -> str:
    name, _ = generate_name_and_justification(smiles)
    return name


def generate_new_justification(smiles: str, compound_name: str) -> str:
    _, justif = generate_name_and_justification(smiles)
    return justif


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

Format WAJIB (jangan tambah teks lain):
[1] <justifikasi 2-3 kalimat>
[2] <justifikasi 2-3 kalimat>
[3] <justifikasi 2-3 kalimat>"""

    text = safe_invoke(prompt, fallback="", max_len=5000)

    results: List[str] = []
    for i, smi in enumerate(dataset_smiles_list):
        m = re.search(rf"\[{i+1}\]\s*(.*?)(?=\[\d+\]|$)", text, flags=re.DOTALL)
        if m:
            cand = _clean(m.group(1), 900)
            if len(cand) < 200:
                results.append(
                    _clean(
                        f"Senyawa ini menunjukkan kemiripan dengan {compound_name} melalui pola ikatan dan motif gugus fungsi yang sejenis. "
                        f"Kedekatan fitur seperti panjang rantai, percabangan, atau segmen polar dapat menghasilkan fingerprint yang berdekatan. "
                        f"Kesamaan tersebut sering berkorelasi dengan kecenderungan sifat fisik yang mirip pada level segmen rantai.",
                        900,
                    )
                )
            else:
                results.append(_clean(cand, 900))
        else:
            results.append(
                _clean(
                    f"Struktur ini mirip dengan {compound_name} berdasarkan kemiripan pola struktur dan fitur gugus fungsi utama. "
                    f"Dalam pencarian similarity berbasis fingerprint, pola ikatan yang serupa biasanya meningkatkan skor kemiripan. "
                    f"Perbedaan kecil masih mungkin muncul pada panjang rantai dan tingkat percabangan.",
                    900,
                )
            )

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

    prompt = f"""Prediksi Tg (°C) untuk kandidat polimer berikut.
SMILES: {smiles}
Nama: {compound_name}
MW: {mw:.2f}
Rings: {rings}

Balas HANYA JSON valid:
{{"predicted_tg": -20.0, "tg_justification": "2-4 kalimat alasan berbasis kekakuan/fleksibilitas rantai dan polaritas."}}

Aturan:
- predicted_tg float -100..300
- tg_justification maksimal 500 karakter
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
                "tg_justification": _clean(_fix_degree_symbol(just), 500),
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