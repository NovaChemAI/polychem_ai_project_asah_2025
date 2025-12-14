import os
import time
import pandas as pd
from typing import Optional, List
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import requests
from io import StringIO
from typing import Optional, List

RDLogger.DisableLog("rdApp.*")

# GLOBAL CACHED STATE
df: Optional[pd.DataFrame] = None
dataset_rdkit_fps: Optional[List] = None  # ExplicitBitVect list

NBITS = 2048
RADIUS = 2

# ID Google Drive
ID_DATA = "1TbnPcrxCysz-eOaSfoLfBxPohBQq7JaY"

# local cache
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_CACHE_DIR = os.path.join(BASE_DIR, "data_cache")
os.makedirs(DATA_CACHE_DIR, exist_ok=True)
LOCAL_CSV_PATH = os.path.join(DATA_CACHE_DIR, "dataset.csv")


def _download_drive_csv_text(file_id: str, timeout: int = 60) -> str:
    """
    Download CSV dari Google Drive (robust):
    - user-agent
    - handle confirm token warning page
    Return: raw CSV text
    """
    url = "https://drive.google.com/uc?export=download"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                      "AppleWebKit/537.36 (KHTML, like Gecko) "
                      "Chrome/120.0 Safari/537.36"
    }

    with requests.Session() as s:
        resp = s.get(url, params={"id": file_id}, headers=headers, timeout=timeout)
        resp.raise_for_status()

        # Kadang Drive kasih HTML warning + cookie confirm
        ctype = (resp.headers.get("Content-Type") or "").lower()
        if "text/html" in ctype:
            confirm_token = None
            for k, v in resp.cookies.items():
                if k.startswith("download_warning"):
                    confirm_token = v
                    break

            if confirm_token:
                resp = s.get(
                    url,
                    params={"id": file_id, "confirm": confirm_token},
                    headers=headers,
                    timeout=timeout,
                )
                resp.raise_for_status()

        return resp.text


def _download_drive_csv(file_id: str, retries: int = 4, timeout: int = 60) -> pd.DataFrame:
    """
    Download CSV dengan retry + backoff.
    Jika berhasil, simpan ke cache lokal.
    """
    last_err: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            text = _download_drive_csv_text(file_id, timeout=timeout)

            # quick sanity: harus ada kolom SMILES di header
            # (tidak 100% wajib, tapi bantu deteksi HTML error page)
            if "SMILES" not in text.splitlines()[0]:
                # kalau header bukan CSV yang bener, biar dianggap gagal
                raise RuntimeError("Response bukan CSV valid (header tidak mengandung 'SMILES').")

            # simpan cache lokal
            with open(LOCAL_CSV_PATH, "w", encoding="utf-8", newline="") as f:
                f.write(text)

            return pd.read_csv(StringIO(text))

        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 12)  # backoff max 12s
            print(f"⚠️ Download dataset gagal (attempt {attempt}/{retries}): {repr(e)}")
            if attempt < retries:
                time.sleep(wait)

    raise RuntimeError(f"Gagal download dataset dari Google Drive setelah {retries}x: {last_err!r}")

def _load_local_csv_if_exists() -> Optional[pd.DataFrame]:
    """
    Load dataset dari cache lokal kalau ada.
    """
    if os.path.exists(LOCAL_CSV_PATH):
        try:
            print(f"📦 Fallback: load dataset dari cache lokal: {LOCAL_CSV_PATH}")
            return pd.read_csv(LOCAL_CSV_PATH)
        except Exception as e:
            print("❌ Gagal baca cache lokal:", repr(e))
            return None
    return None


def load_dataset() -> None:
    """
    Load dataset saat startup:
    1) Coba download dari Drive (retry)
    2) Kalau gagal, fallback ke cache lokal
    3) Precompute RDKit fingerprints
    """
    global df, dataset_rdkit_fps

    print("🔥 load_dataset() jalan...")

    raw = None
    try:
        raw = _download_drive_csv(ID_DATA, retries=4, timeout=60)
        print("✅ Download dataset dari Drive: OK")
    except Exception as e:
        print("⚠️ Download Drive gagal, coba pakai cache lokal...")
        raw = _load_local_csv_if_exists()
        if raw is None:
            # pilihan: kalau benar-benar tidak ada dataset, tetap jalan tapi kosong
            print("❌ Tidak ada dataset (Drive gagal & cache lokal tidak ada). Server tetap jalan, tapi similarity kosong.")
            df = pd.DataFrame(columns=["SMILES"])
            dataset_rdkit_fps = []
            return

    if "SMILES" not in raw.columns:
        raise RuntimeError("Dataset tidak punya kolom 'SMILES'.")

    print(f"📌 Columns: {list(raw.columns)}")

    valid_rows = []
    fps_rdkit: List = []

    for row in raw.itertuples(index=False):
        smi = getattr(row, "SMILES", None)
        if not isinstance(smi, str) or not smi.strip():
            continue

        smi = smi.strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=NBITS)
        fps_rdkit.append(fp)
        valid_rows.append(row._asdict())

    df = pd.DataFrame(valid_rows)
    dataset_rdkit_fps = fps_rdkit

    print(f"✅ Dataset loaded: {len(df)} baris valid")
    print(f"✅ Fingerprints ready (rdkit): {len(dataset_rdkit_fps)}")