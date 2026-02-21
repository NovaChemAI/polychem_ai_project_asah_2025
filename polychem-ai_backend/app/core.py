import numpy as np
from functools import lru_cache
from typing import List, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

import app.store as store
from app.llm import (
    generate_compound_name,
    generate_new_justification,
    justify_similar_compounds_batch,
    predict_tg_with_llm,
)

# OPTIONAL: kalau fungsi ini ada di llm.py, kita pakai.
# Kalau tidak ada, core.py tetap jalan (fallback internal).
try:
    from app.llm import tg_fallback_heuristic as _tg_fallback_heuristic
except Exception:
    _tg_fallback_heuristic = None

from app.images import save_smiles_png
from app.cache import cache, key_new_compound, key_similar_justif, push_history

from app.settings import COMPOUNDS_DIR  # ✅ Koyeb-safe (default: /tmp/static/compounds)


# ============================================================
# UTILS
# ============================================================

def build_image_url(filename: str) -> str:
    return f"/static/compounds/{filename}"


def normalize_smiles(smiles: str) -> str:
    return (smiles or "").strip()


def is_too_simple(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return True
    return mol.GetNumAtoms() <= 5


def build_fingerprints(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    fp_rdkit = AllChem.GetMorganFingerprintAsBitVect(
        mol, store.RADIUS, nBits=store.NBITS
    )
    fp_list = list(fp_rdkit)
    return fp_list, fp_rdkit


def _tg_fallback_local(smiles: str) -> float:
    """
    Heuristik cepat berbasis RDKit supaya Tg similar compounds tidak jadi 0.
    (Tidak akurat ilmiah, tapi lebih masuk akal daripada 0 terus.)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0

    mw = Descriptors.MolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)

    tg = 10.0
    tg += 18.0 * rings
    tg -= 6.0 * rot
    tg += 8.0 * (hbd + 0.6 * hba)

    if mw > 200:
        tg += 20.0
    elif mw > 100:
        tg += 10.0

    tg = max(-80.0, min(250.0, tg))
    return round(float(tg), 1)


def _get_metadata_for_smiles(smiles: str) -> dict:
    """
    Ambil metadata dari dataset untuk SMILES tertentu.
    Kalau tidak ada → fallback aman.
    """
    fallback = {
        "name": "",
        "formula": "",
        "molecular_weight": 0.0,
        "tg": 0.0,
        "pid": "",
        "polymer_class": "",
    }

    if store.df is None or "SMILES" not in store.df.columns:
        return fallback

    s = str(smiles or "").strip()
    row = store.df[store.df["SMILES"].astype(str).str.strip() == s]
    if row.empty:
        return fallback

    r = row.iloc[0]

    def pick(keys, default):
        for k in keys:
            if k in row.columns:
                v = r.get(k, default)
                if v is not None and v != "":
                    return v
        return default

    name = str(pick(["Name", "IUPAC_Name"], "") or "")
    formula = str(pick(["Formula", "Molecular_Formula"], "") or "")

    mw = pick(["MolecularWeight", "MolecularWeight(g/mol)", "Molecular_Weight"], 0.0) or 0.0
    tg = pick(["Tg", "TG", "Tg (°C)"], 0.0) or 0.0

    pid = str(pick(["PID"], "") or "")
    polymer_class = str(pick(["Polymer Class", "PolymerClass", "Polymer_Class"], "") or "")

    try:
        mw = float(mw)
    except Exception:
        mw = 0.0

    try:
        tg = float(tg)
    except Exception:
        tg = 0.0

    return {
        "name": name,
        "formula": formula,
        "molecular_weight": mw,
        "tg": tg,
        "pid": pid,
        "polymer_class": polymer_class,
    }


# ============================================================
# CACHING (L1 RAM + L2 diskcache)
# ============================================================

@lru_cache(maxsize=256)
def _cached_new_compound_llm(smiles_norm: str):
    dkey = key_new_compound(smiles_norm)
    cached = cache.get(dkey)
    if isinstance(cached, dict) and "name" in cached and "justifikasi" in cached:
        return cached["name"], cached["justifikasi"]

    name = generate_compound_name(smiles_norm)
    justif = generate_new_justification(smiles_norm, name)

    cache.set(dkey, {"name": name, "justifikasi": justif}, expire=60 * 60 * 24 * 7)
    return name, justif


@lru_cache(maxsize=256)
def _cached_similar_justifs(compound_name: str, top_smiles: Tuple[str, ...]) -> List[str]:
    dkey = key_similar_justif(compound_name, list(top_smiles))
    cached = cache.get(dkey)
    if isinstance(cached, list):
        return cached

    justifs = justify_similar_compounds_batch(compound_name, list(top_smiles))
    if not isinstance(justifs, list):
        justifs = ["Mirip secara struktur umum (fallback)."] * len(top_smiles)

    cache.set(dkey, justifs, expire=60 * 60 * 24)
    return justifs


# ============================================================
# SIMILARITY SEARCH
# ============================================================

def find_similar_compounds(input_fp_rdkit, compound_name: str, top_k: int = 3):
    if store.df is None or store.dataset_rdkit_fps is None:
        return []

    similarities = DataStructs.BulkTanimotoSimilarity(
        input_fp_rdkit, store.dataset_rdkit_fps
    )
    similarities = np.array(similarities, dtype=float)

    top_indices = np.argsort(similarities)[::-1][:top_k]
    top_smiles = tuple(str(store.df.iloc[idx]["SMILES"]).strip() for idx in top_indices)

    justifs = _cached_similar_justifs(compound_name, top_smiles)
    if not isinstance(justifs, list):
        justifs = ["Mirip secara struktur umum (fallback)."] * len(top_smiles)

    similar_compounds = []
    for i, idx in enumerate(top_indices):
        ds_smiles = top_smiles[i]
        meta = _get_metadata_for_smiles(ds_smiles)

        mol_ds = Chem.MolFromSmiles(ds_smiles)
        if mol_ds is not None:
            if not meta["formula"]:
                meta["formula"] = rdMolDescriptors.CalcMolFormula(mol_ds)

            if meta["molecular_weight"] <= 0.0:
                meta["molecular_weight"] = round(float(Descriptors.MolWt(mol_ds)), 2)

        # ✅ pastikan Tg untuk rank 1/2/3 tidak 0
        if float(meta.get("tg", 0.0)) == 0.0:
            if _tg_fallback_heuristic is not None:
                try:
                    est = _tg_fallback_heuristic(ds_smiles)
                    meta["tg"] = float(est.get("tg", 0.0)) if isinstance(est, dict) else float(est)
                except Exception:
                    meta["tg"] = _tg_fallback_local(ds_smiles)
            else:
                meta["tg"] = _tg_fallback_local(ds_smiles)

        filename = ""
        image_url = ""
        try:
            filename = save_smiles_png(ds_smiles, COMPOUNDS_DIR)
            image_url = build_image_url(filename)
        except Exception:
            filename = ""
            image_url = ""
        score = float(similarities[idx])

        # Fallback nama: pakai formula atau polymer_class kalau name kosong
        display_name = meta["name"]
        if not display_name or display_name in ["DummyName", "GeneratedCompound"]:
            if meta["polymer_class"]:
                display_name = f"{meta['polymer_class']} ({meta['formula']})"
            elif meta["formula"]:
                display_name = f"Polymer {meta['formula']}"
            else:
                display_name = "Unknown Polymer"

        similar_compounds.append({
            "rank": i + 1,
            "smiles": ds_smiles,
            "name": display_name,
            "formula": meta["formula"],
            "molecular_weight": meta["molecular_weight"],
            "tg": float(meta["tg"]),
            "pid": meta["pid"],
            "polymer_class": meta["polymer_class"],
            "similarity_score": score,
            "similarity_percent": score * 100.0,
            "justifikasi": justifs[i] if i < len(justifs) else "Mirip secara struktur umum (fallback).",
            "image_filename": filename,
            "image_url": image_url,
        })

    return similar_compounds


# ============================================================
# ORCHESTRATOR
# ============================================================

def recommend_new_compound(input_smiles: str) -> dict:
    smiles_norm = normalize_smiles(input_smiles)

    fp_list, input_fp_rdkit = build_fingerprints(smiles_norm)
    if input_fp_rdkit is None:
        return {"error": "SMILES tidak valid", "input": input_smiles}

    meta_input = _get_metadata_for_smiles(smiles_norm)

    # isi formula & mw dari RDKit jika kosong
    mol = Chem.MolFromSmiles(smiles_norm)
    if mol is not None:
        if not meta_input["formula"]:
            meta_input["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        if meta_input["molecular_weight"] <= 0.0:
            meta_input["molecular_weight"] = round(float(Descriptors.MolWt(mol)), 2)

    # name + justifikasi
    tg_justification = ""
    if is_too_simple(smiles_norm):
        compound_name = meta_input["name"] or "Simple Compound"
        justifikasi = "SMILES ini sangat sederhana dan kemungkinan besar sudah umum."
        tg_justification = "Tg tidak diprediksi untuk senyawa sangat sederhana."
    else:
        if meta_input["name"] and meta_input["name"] not in ["DummyName", "GeneratedCompound"]:
            compound_name = meta_input["name"]
            justifikasi = generate_new_justification(smiles_norm, compound_name)
        else:
            compound_name, justifikasi = _cached_new_compound_llm(smiles_norm)

            if compound_name in ["DummyName", "GeneratedCompound", "Unknown", ""]:
                # Fallback: pakai polymer_class atau formula
                if meta_input["polymer_class"]:
                    compound_name = f"{meta_input['polymer_class']} ({meta_input['formula']})"
                elif meta_input["formula"]:
                    compound_name = f"Polymer {meta_input['formula']}"
                else:
                    compound_name = "Novel Polymer"

        # Tg input: pakai dataset kalau ada, kalau tidak -> LLM
        if float(meta_input.get("tg", 0.0)) == 0.0:
            tg_data = predict_tg_with_llm(smiles_norm, compound_name)
            meta_input["tg"] = float(tg_data.get("tg", 0.0))
            tg_justification = str(tg_data.get("tg_justification", ""))
        else:
            tg_justification = "Data Tg dari database (exact match)."

    new_filename = ""
    new_image_url = ""
    try:
        new_filename = save_smiles_png(smiles_norm, COMPOUNDS_DIR)
        new_image_url = build_image_url(new_filename)
    except Exception:
        new_filename = ""
        new_image_url = ""

    similar_compounds = find_similar_compounds(
        input_fp_rdkit,
        compound_name,
        top_k=3
    )

    result = {
        "status": "success",
        "input_smiles": smiles_norm,
        "new_compound": {
            "name": compound_name,
            "smiles": smiles_norm,
            "formula": meta_input["formula"],
            "molecular_weight": meta_input["molecular_weight"],
            "tg_justification": tg_justification,
            "tg": float(meta_input["tg"]),
            "pid": meta_input["pid"],
            "polymer_class": meta_input["polymer_class"],
            "justifikasi": justifikasi,
            "fingerprint_length": len(fp_list) if fp_list else 0,
            "image_filename": new_filename,
            "image_url": new_image_url,
        },
        "similar_compounds": similar_compounds,
    }

    push_history({"input_smiles": smiles_norm, "result": result})
    return result