from app.firebase_admin_config import db
from firebase_admin import firestore
import hashlib

COLLECTION_NAME = "saved_chemicals"

def _smiles_to_id(smiles: str) -> str:
    return hashlib.md5(smiles.encode()).hexdigest()[:16]

def save_to_library(uid: str, data: dict) -> bool:
    smiles = data.get("smiles", "")
    if data.get("id") and not data.get("isAiResult"):
        doc_id = str(data.get("id"))
    else:
        doc_id = f"ai_{_smiles_to_id(smiles)}"

    properties = data.get("properties", {})
    if isinstance(properties, dict):
        properties = str(properties)

    image_url = data.get("image_url") or data.get("image") or ""

    payload = {
        "id": doc_id,
        "userId": uid,
        "name": data.get("name", "Unknown Compound"),
        "smiles": smiles,
        "category": data.get("category", "Uncategorized"),
        "image": image_url,
        "image_url": image_url,
        "properties": properties,
        "score": data.get("score", "N/A"),
        "isAiResult": data.get("isAiResult", True),
        "savedAt": firestore.SERVER_TIMESTAMP,
    }

    composite_id = f"{uid}_{doc_id}"
    db.collection(COLLECTION_NAME).document(composite_id).set(payload)
    return True


def get_user_library(uid: str) -> list:
    docs = db.collection(COLLECTION_NAME).where("userId", "==", uid).stream()
    return [doc.to_dict() for doc in docs]


def remove_from_library(uid: str, item_id: str) -> bool:
    composite_id = f"{uid}_{item_id}"
    db.collection(COLLECTION_NAME).document(composite_id).delete()
    return True


def check_is_saved_by_smiles(uid: str, smiles: str) -> dict:
    doc_id = f"ai_{_smiles_to_id(smiles)}"
    composite_id = f"{uid}_{doc_id}"
    doc = db.collection(COLLECTION_NAME).document(composite_id).get()
    return {"isSaved": doc.exists, "itemId": doc_id}