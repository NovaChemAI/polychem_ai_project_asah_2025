from app.firebase_admin_config import db
from firebase_admin import firestore

USERS_COLLECTION = "users"

def sync_user_profile(uid: str, name: str, email: str, photo_url: str = "") -> dict:
    user_ref = db.collection(USERS_COLLECTION).document(uid)
    user_doc = user_ref.get()

    if not user_doc.exists:
        payload = {
            "uid": uid,
            "name": name,
            "email": email,
            "role": "user",
            "createdAt": firestore.SERVER_TIMESTAMP,
            "photoURL": photo_url,
        }
        user_ref.set(payload)

        # Ambil ulang dari Firestore supaya createdAt sudah jadi timestamp asli, bukan sentinel
        return user_ref.get().to_dict()

    return user_doc.to_dict()