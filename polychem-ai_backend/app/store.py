import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import requests
from io import StringIO
from typing import Optional, List

RDLogger.DisableLog("rdApp.*")

# GLOBAL CACHED STATE
df: Optional[pd.DataFrame] = None
dataset_rdkit_fps: Optional[List] = None   # ExplicitBitVect list 

NBITS = 2048
RADIUS = 2

# ID Google Drive
ID_DATA = "1TbnPcrxCysz-eOaSfoLfBxPohBQq7JaY"


def _download_drive_csv(file_id: str) -> pd.DataFrame:
    """
    Download CSV dari Google Drive dengan cara stabil (handle warning page).
    Return DataFrame.
    """
    url = "https://drive.google.com/uc?export=download"

    try:
        with requests.Session() as s:
            resp = s.get(url, params={"id": file_id}, timeout=60)
            resp.raise_for_status()
            if "text/html" in resp.headers.get("Content-Type", ""):
                for k, v in resp.cookies.items():
                    if k.startswith("download_warning"):
                        resp = s.get(
                            url, params={"id": file_id, "confirm": v}, timeout=60
                        )
                        resp.raise_for_status()
                        break

            return pd.read_csv(StringIO(resp.text))

    except requests.RequestException as e:
        raise RuntimeError(f"Gagal download dataset dari Google Drive: {e}")


def load_dataset() -> None:
    """
    Load dataset dari Google Drive lalu precompute fingerprint sekali saat startup.

    Workflow:
    1) Download CSV (robust)
    2) Pastikan ada kolom SMILES
    3) Filter SMILES invalid
    4) Hitung Morgan fingerprint (RDKit ExplicitBitVect)
    5) Cache df + dataset_rdkit_fps
    """
    global df, dataset_rdkit_fps

    print("🔥 load_dataset() jalan...")

    raw = _download_drive_csv(ID_DATA)

    if "SMILES" not in raw.columns:
        raise RuntimeError("Dataset tidak punya kolom 'SMILES'.")

    valid_rows = []
    fps_rdkit: List = []

    # itertuples lebih cepat dari iterrows
    for row in raw.itertuples(index=False):
        smi = getattr(row, "SMILES", None)

        if not isinstance(smi, str) or not smi.strip():
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=NBITS
        )

        fps_rdkit.append(fp)
        valid_rows.append(row._asdict())

    df = pd.DataFrame(valid_rows)
    dataset_rdkit_fps = fps_rdkit

    print(f"✅ Dataset loaded: {len(df)} baris valid")
    print(f"✅ Fingerprints ready (rdkit): {len(dataset_rdkit_fps)}")
