from app.firebase_admin_config import db
from firebase_admin import firestore
import json

HISTORY_COLLECTION = "search_history"

def add_to_history(uid: str, query_text: str, result_data: dict):
    nc = result_data.get("new_compound", {})
    payload = {
        "userId": uid,
        "query": query_text,
        "result_name": nc.get("name", "Unknown"),
        "result_smiles": nc.get("smiles", ""),
        "full_data": json.dumps(result_data),
        "timestamp": firestore.SERVER_TIMESTAMP,
    }
    db.collection(HISTORY_COLLECTION).add(payload)


def get_user_history(uid: str, limit_count: int = 50) -> list:
    docs = (
        db.collection(HISTORY_COLLECTION)
        .where("userId", "==", uid)
        .order_by("timestamp", direction=firestore.Query.DESCENDING)
        .limit(limit_count)
        .stream()
    )
    items = []
    for doc in docs:
        d = doc.to_dict()
        d["id"] = doc.id
        items.append(d)
    return items