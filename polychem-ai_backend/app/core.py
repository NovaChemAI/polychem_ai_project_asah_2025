import os
import numpy as np
from functools import lru_cache
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from typing import List, Tuple

import app.store as store
from app.llm import (
    generate_compound_name,
    generate_new_justification,
    justify_similar_compounds_batch
)
from app.images import save_smiles_png

# disk cache
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
    """
    Ambil metadata dari dataset untuk SMILES tertentu.
    Kalau tidak ada di dataset → fallback aman.
    """
    fallback = {
        "name": "",
        "formula": "",
        "molecular_weight": 0.0,
        "tg": 0.0,
        "pid": "",
        "polymer_class": "",
    }

    if store.df is None:
        return fallback

    # cari exact match
    row = store.df[store.df["SMILES"] == smiles]
    if row.empty:
        return fallback

    r = row.iloc[0]
    return {
        "name": str(r.get("Name", "")),
        "formula": str(r.get("Formula", "")),
        "molecular_weight": float(r.get("MolecularWeight", 0.0) or 0.0),
        "tg": float(r.get("Tg", 0.0) or 0.0),
        "pid": str(r.get("PID", "")),
        "polymer_class": str(r.get("Polymer Class", "")),
    }


# ============================================================
# L1 cache (RAM) + L2 cache (diskcache)
# ============================================================

@lru_cache(maxsize=256)
def _cached_new_compound_llm(smiles_norm: str):
    dkey = key_new_compound(smiles_norm)
    cached = cache.get(dkey)
    if cached is not None:
        return cached["name"], cached["justifikasi"]

    name = generate_compound_name(smiles_norm)
    justif = generate_new_justification(smiles_norm, name)

    cache.set(
        dkey,
        {"name": name, "justifikasi": justif},
        expire=60 * 60 * 24 * 7  # 7 hari
    )
    return name, justif


@lru_cache(maxsize=256)
def _cached_similar_justifs(compound_name: str, top_smiles: Tuple[str, ...]) -> List[str]:
    dkey = key_similar_justif(compound_name, list(top_smiles))
    cached = cache.get(dkey)
    if cached is not None:
        return cached

    justifs = justify_similar_compounds_batch(compound_name, list(top_smiles))
    cache.set(
        dkey,
        justifs,
        expire=60 * 60 * 24  # 1 hari
    )
    return justifs


# ============================================================
# Similarity search
# ============================================================

def find_similar_compounds(input_fp_rdkit, compound_name: str, top_k: int = 3):
    if store.df is None or store.dataset_rdkit_fps is None:
        return []

    similarities = DataStructs.BulkTanimotoSimilarity(
        input_fp_rdkit, store.dataset_rdkit_fps
    )
    similarities = np.array(similarities, dtype=float)

    top_indices = np.argsort(similarities)[::-1][:top_k]
    top_smiles = tuple(store.df.iloc[idx]["SMILES"] for idx in top_indices)

    justifs = _cached_similar_justifs(compound_name, top_smiles)

    similar_compounds = []
    for i, idx in enumerate(top_indices):
        ds_smiles = top_smiles[i]
        meta = _get_metadata_for_smiles(ds_smiles)

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
            "justifikasi": justifs[i] if i < len(justifs) else "",
            "image_filename": filename,
            "image_url": build_image_url(filename),
        })

    return similar_compounds


# ============================================================
# Orchestrator
# ============================================================

def recommend_new_compound(input_smiles: str) -> dict:
    smiles_norm = normalize_smiles(input_smiles)

    fp_list, input_fp_rdkit = build_fingerprints(smiles_norm)
    if input_fp_rdkit is None:
        return {"error": "SMILES tidak valid", "input": input_smiles}

    # metadata dari dataset jika ada
    meta_input = _get_metadata_for_smiles(smiles_norm)

    # kalau terlalu sederhana, skip LLM
    if is_too_simple(smiles_norm):
        compound_name = meta_input["name"] or "Known simple compound"
        justifikasi = (
            "SMILES ini sangat sederhana dan kemungkinan besar adalah senyawa umum "
            "yang sudah teridentifikasi luas, jadi tidak dibuat nama baru."
        )
    else:
        # kalau dataset punya name, pakai itu dulu (lebih masuk akal)
        if meta_input["name"]:
            compound_name = meta_input["name"]
            justifikasi = generate_new_justification(smiles_norm, compound_name)
        else:
            compound_name, justifikasi = _cached_new_compound_llm(smiles_norm)

    new_filename = save_smiles_png(smiles_norm, COMPOUND_DIR)

    similar_compounds = find_similar_compounds(
        input_fp_rdkit,
        compound_name,
        top_k=3
    )

    # ✅ hasil final (urut field dijaga)
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
            "justifikasi": justifikasi,
            "fingerprint_length": len(fp_list) if fp_list else 0,
            "image_filename": new_filename,
            "image_url": build_image_url(new_filename),
        },
        "similar_compounds": similar_compounds
    }

    # ✅ simpan history (diskcache)
    push_history({
        "input_smiles": smiles_norm,
        "result": result
    })

    return result
