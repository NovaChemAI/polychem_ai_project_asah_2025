import os
import hashlib
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image as PILImage


def smiles_to_image(smiles: str, image_size=(800, 300)) -> Optional[PILImage.Image]:
    """SMILES -> PIL Image. Return None kalau invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    AllChem.Compute2DCoords(mol)

    return Draw.MolToImage(
        mol,
        size=image_size,
        kekulize=True,
        wedgeBonds=True,
    )


def get_compounds_dir() -> str:
    """
    Path absolut ke static/compounds (root project).
    app/images.py -> naik 1 level ke root -> static/compounds
    """
    base_dir = os.path.dirname(os.path.dirname(__file__))
    return os.path.join(base_dir, "static", "compounds")


def save_smiles_png(smiles: str, out_dir: Optional[str] = None) -> str:
    """
    Simpan PNG ke out_dir (default static/compounds).
    Cache by hash: kalau file sudah ada, gak generate ulang.
    Return filename saja.
    """
    smiles = smiles.strip()

    if out_dir is None:
        out_dir = get_compounds_dir()

    os.makedirs(out_dir, exist_ok=True)

    hash8 = hashlib.md5(smiles.encode("utf-8")).hexdigest()[:8]
    filename = f"compound_{hash8}.png"
    file_path = os.path.join(out_dir, filename)

    if os.path.exists(file_path):
        return filename  # cache hit

    img = smiles_to_image(smiles)
    if img is None:
        raise ValueError("SMILES tidak valid, gambar gagal dibuat")

    img.save(file_path, format="PNG")
    return filename
