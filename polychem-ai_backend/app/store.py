import os
import time
from io import StringIO
from typing import Optional, List

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# =============================================================================
# GLOBAL CACHED STATE
# =============================================================================
df: Optional[pd.DataFrame] = None
dataset_rdkit_fps: Optional[List] = None  # list of ExplicitBitVect

NBITS = 2048
RADIUS = 2

# ID Google Drive
ID_DATA = os.getenv("ID_DATASET_DRIVE", "13L5ZFx_vZyrTwS4tUeWADO1-SE_oPHkM")

# =============================================================================
# LOCAL CACHE (Koyeb-safe)
# Default /tmp writable. Bisa override pakai env DATA_CACHE_DIR
# =============================================================================
DATA_CACHE_DIR = os.getenv("DATA_CACHE_DIR", "/tmp/data_cache")
os.makedirs(DATA_CACHE_DIR, exist_ok=True)
LOCAL_CSV_PATH = os.path.join(DATA_CACHE_DIR, "dataset.csv")


def _download_drive_csv_text(file_id: str, timeout: int = 60) -> str:
    """
    Download CSV dari Google Drive (robust):
    - User-Agent
    - handle confirm token warning page
    Return: raw CSV text
    """
    url = "https://drive.google.com/uc?export=download"
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0 Safari/537.36"
        )
    }

    with requests.Session() as s:
        resp = s.get(url, params={"id": file_id}, headers=headers, timeout=timeout)
        resp.raise_for_status()

        # Drive kadang ngasih HTML warning + cookie confirm
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


def _looks_like_csv(text: str) -> bool:
    """
    Heuristik cepat biar gak salah parse HTML error page.
    """
    if not text:
        return False
    head = "\n".join(text.splitlines()[:5]).lower()
    if "<html" in head or "<!doctype html" in head:
        return False
    # minimal ada comma/semicolon dan ada baris pertama
    first = (text.splitlines()[0] if text.splitlines() else "").strip()
    if not first:
        return False
    if ("," not in first) and (";" not in first) and ("\t" not in first):
        # masih mungkin CSV 1 kolom, tapi dataset kamu harusnya multi kolom
        return False
    return True


def _save_local_cache(text: str) -> None:
    """
    Simpan dataset ke cache lokal. Kalau gagal (FS read-only), jangan crash.
    """
    try:
        with open(LOCAL_CSV_PATH, "w", encoding="utf-8", newline="") as f:
            f.write(text)
    except Exception as e:
        print("⚠️ Gagal simpan cache dataset ke disk:", repr(e))


def _download_drive_csv(file_id: str, retries: int = 4, timeout: int = 60) -> pd.DataFrame:
    """
    Download CSV dengan retry + backoff.
    Jika berhasil, simpan ke cache lokal (/tmp).
    """
    last_err: Optional[Exception] = None

    for attempt in range(1, retries + 1):
        try:
            text = _download_drive_csv_text(file_id, timeout=timeout)

            if not _looks_like_csv(text):
                raise RuntimeError("Response bukan CSV valid (kemungkinan halaman HTML / error).")

            # parse CSV
            raw = pd.read_csv(StringIO(text))

            # simpan cache lokal (best effort)
            _save_local_cache(text)

            return raw

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


def _ensure_smiles_column(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Pastikan ada kolom 'SMILES' persis.
    Kadang dataset punya 'smiles' atau 'Smiles' dsb.
    """
    if raw is None or raw.empty:
        return pd.DataFrame(columns=["SMILES"])

    cols = list(raw.columns)
    if "SMILES" in cols:
        return raw

    # coba cari kolom yang mirip
    lowered = {c.lower().strip(): c for c in cols}
    for key in ["smiles", "smile"]:
        if key in lowered:
            raw = raw.rename(columns={lowered[key]: "SMILES"})
            return raw

    return raw  # biar error ditangani di load_dataset()


def load_dataset() -> None:
    """
    Load dataset saat startup:
    1) Coba download dari Drive (retry)
    2) Kalau gagal, fallback ke cache lokal (/tmp)
    3) Precompute RDKit fingerprints
    """
    global df, dataset_rdkit_fps

    print("🔥 load_dataset() jalan...")

    raw: Optional[pd.DataFrame] = None

    try:
        raw = _download_drive_csv(ID_DATA, retries=4, timeout=60)
        print("✅ Download dataset dari Drive: OK")
    except Exception as e:
        print("⚠️ Download Drive gagal, coba pakai cache lokal...")
        print("   sebab:", repr(e))
        raw = _load_local_csv_if_exists()

        if raw is None:
            print("❌ Tidak ada dataset (Drive gagal & cache lokal tidak ada).")
            print("   Server tetap jalan, tapi similarity kosong.")
            df = pd.DataFrame(columns=["SMILES"])
            dataset_rdkit_fps = []
            return

    raw = _ensure_smiles_column(raw)

    if "SMILES" not in raw.columns:
        # jangan bikin server mati total; fallback kosong biar /predict tetap jalan
        print("❌ Dataset tidak punya kolom 'SMILES'. Server tetap jalan, similarity kosong.")
        df = pd.DataFrame(columns=["SMILES"])
        dataset_rdkit_fps = []
        return

    print(f"📌 Columns: {list(raw.columns)}")

    valid_rows: List[dict] = []
    fps_rdkit: List = []

    # iterasi per row; aman untuk dataset besar
    for row in raw.to_dict(orient="records"):
        smi = row.get("SMILES", None)
        if not isinstance(smi, str) or not smi.strip():
            continue

        smi = smi.strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=NBITS)
        fps_rdkit.append(fp)
        valid_rows.append(row)

    df = pd.DataFrame(valid_rows)
    dataset_rdkit_fps = fps_rdkit

    print(f"✅ Dataset loaded: {len(df)} baris valid")
    print(f"✅ Fingerprints ready (rdkit): {len(dataset_rdkit_fps)}")