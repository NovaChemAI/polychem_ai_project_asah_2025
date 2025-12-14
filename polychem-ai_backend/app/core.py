import os
import numpy as np
from functools import lru_cache
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from typing import List, Tuple

import app.store as store
from app.llm import (
    generate_compound_name,
    generate_new_justification,
    justify_similar_compounds_batch,
    predict_tg_with_llm,
    tg_fallback_heuristic, # <--- IMPORT BERHASIL JIKA LLM.PY SUDAH DIPERBAIKI
)
from app.images import save_smiles_png
from app.cache import cache, key_new_compound, key_similar_justif, push_history

# path absolut ke static/compounds
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
COMPOUND_DIR = os.path.join(BASE_DIR, "static", "compounds")


def build_image_url(filename: str) -> str:
    return f"/static/compounds/{filename}"


def normalize_smiles(smiles: str) -> str:
    return smiles.strip()


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


def _get_metadata_for_smiles(smiles: str) -> dict:
    fallback = {
        "name": "", "formula": "", "molecular_weight": 0.0,
        "tg": 0.0, "pid": "", "polymer_class": "",
    }
    if store.df is None: return fallback
    
    row = store.df[store.df["SMILES"].astype(str).str.strip() == smiles.strip()]
    if row.empty: return fallback

    r = row.iloc[0]
    
    def get_val(keys, default):
        for k in keys:
            if k in r and r[k]: return r[k]
        return default

    mw = get_val(["MolecularWeight", "MolecularWeight(g/mol)"], 0.0)
    tg = get_val(["Tg", "TG", "Tg (°C)"], 0.0)
    
    try: mw = float(mw)
    except: mw = 0.0
    try: tg = float(tg)
    except: tg = 0.0

    return {
        "name": str(get_val(["Name", "IUPAC_Name"], "")),
        "formula": str(get_val(["Formula", "Molecular_Formula"], "")),
        "molecular_weight": mw,
        "tg": tg,
        "pid": str(get_val(["PID"], "")),
        "polymer_class": str(get_val(["Polymer Class"], "")),
    }


# ============================================================
# CACHING
# ============================================================
@lru_cache(maxsize=256)
def _cached_new_compound_llm(smiles_norm: str):
    dkey = key_new_compound(smiles_norm)
    cached = cache.get(dkey)
    if cached: return cached["name"], cached["justifikasi"]

    name = generate_compound_name(smiles_norm)
    justif = generate_new_justification(smiles_norm, name)
    cache.set(dkey, {"name": name, "justifikasi": justif}, expire=60*60*24*7)
    return name, justif


@lru_cache(maxsize=256)
def _cached_similar_justifs(compound_name: str, top_smiles: Tuple[str, ...]) -> List[str]:
    dkey = key_similar_justif(compound_name, list(top_smiles))
    cached = cache.get(dkey)
    if cached: return cached

    justifs = justify_similar_compounds_batch(compound_name, list(top_smiles))
    cache.set(dkey, justifs, expire=60*60*24)
    return justifs


# ============================================================
# SIMILARITY SEARCH
# ============================================================
def find_similar_compounds(input_fp_rdkit, compound_name: str, top_k: int = 3):
    if store.df is None or store.dataset_rdkit_fps is None: return []

    similarities = DataStructs.BulkTanimotoSimilarity(input_fp_rdkit, store.dataset_rdkit_fps)
    similarities = np.array(similarities, dtype=float)
    top_indices = np.argsort(similarities)[::-1][:top_k]
    top_smiles = tuple(store.df.iloc[idx]["SMILES"] for idx in top_indices)

    justifs = _cached_similar_justifs(compound_name, top_smiles)
    if not isinstance(justifs, list): justifs = ["Mirip."] * len(top_smiles)

    similar_compounds = []
    for i, idx in enumerate(top_indices):
        ds_smiles = str(top_smiles[i]).strip()
        meta = _get_metadata_for_smiles(ds_smiles)
        
        # FIX 1: MW & Formula
        mol_ds = Chem.MolFromSmiles(ds_smiles)
        if mol_ds:
            if not meta["formula"]:
                meta["formula"] = rdMolDescriptors.CalcMolFormula(mol_ds)
            if meta["molecular_weight"] < 1.0:
                meta["molecular_weight"] = Descriptors.MolWt(mol_ds)
            meta["molecular_weight"] = round(meta["molecular_weight"], 2)

        # FIX 2: Tg Similar Compounds (Pakai Heuristic agar cepat)
        if meta["tg"] == 0.0:
            est = tg_fallback_heuristic(ds_smiles)
            meta["tg"] = est["tg"]

        filename = save_smiles_png(ds_smiles, COMPOUND_DIR)
        score = float(similarities[idx])
        
        similar_compounds.append({
            "rank": i + 1,
            "smiles": ds_smiles,
            "name": meta["name"],
            "formula": meta["formula"],
            "molecular_weight": meta["molecular_weight"],
            "tg": meta["tg"],
            "pid": meta["pid"],
            "polymer_class": meta["polymer_class"],
            "similarity_score": score,
            "similarity_percent": score * 100.0,
            "justifikasi": justifs[i] if i < len(justifs) else "Mirip.",
            "image_filename": filename,
            "image_url": build_image_url(filename),
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
    
    # 1. Hitung Fisika Kimia
    mol = Chem.MolFromSmiles(smiles_norm)
    if mol:
        if not meta_input["formula"]:
            meta_input["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        if meta_input["molecular_weight"] < 1.0:
            meta_input["molecular_weight"] = Descriptors.MolWt(mol)
        meta_input["molecular_weight"] = round(meta_input["molecular_weight"], 2)

    # 2. Nama & Anti-DummyName
    tg_justification = ""
    compound_name = ""
    justifikasi = ""

    if is_too_simple(smiles_norm):
        compound_name = meta_input["name"] or "Simple Compound"
        justifikasi = "Molekul sederhana."
        tg_justification = "N/A"
    else:
        if meta_input["name"] and meta_input["name"] != "DummyName":
            compound_name = meta_input["name"]
            justifikasi = generate_new_justification(smiles_norm, compound_name)
        else:
            compound_name, justifikasi = _cached_new_compound_llm(smiles_norm)
            # Guard jika nama masih Dummy/Generated
            blacklist = ["DummyName", "GeneratedCompound", "Unknown", ""]
            if compound_name in blacklist:
                compound_name = f"Polymer ({meta_input['formula']})"

        # 3. Prediksi Tg
        if meta_input["tg"] == 0.0:
            print(f"🔮 Memprediksi Tg via LLM untuk {compound_name}...")
            tg_data = predict_tg_with_llm(smiles_norm, compound_name)
            meta_input["tg"] = tg_data.get("tg", 0.0)
            tg_justification = tg_data.get("tg_justification", "")
        else:
            tg_justification = "Data Tg dari database."

    new_filename = save_smiles_png(smiles_norm, COMPOUND_DIR)
    
    similar_compounds = find_similar_compounds(
        input_fp_rdkit, compound_name, top_k=3
    )

    result = {
        "status": "success",
        "input_smiles": smiles_norm,
        "new_compound": {
            "name": compound_name,
            "smiles": smiles_norm,
            "formula": meta_input["formula"],
            "molecular_weight": meta_input["molecular_weight"],
            "tg": meta_input["tg"],
            "pid": meta_input["pid"],
            "polymer_class": meta_input["polymer_class"],
            "tg_justification": tg_justification,
            "justifikasi": justifikasi,
            "fingerprint_length": len(fp_list) if fp_list else 0,
            "image_filename": new_filename,
            "image_url": build_image_url(new_filename),
        },
        "similar_compounds": similar_compounds
    }

    push_history({"input_smiles": smiles_norm, "result": result})
    return result