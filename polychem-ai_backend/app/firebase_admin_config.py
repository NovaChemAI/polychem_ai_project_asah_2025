import firebase_admin
from firebase_admin import credentials, firestore
import os
import json

# Coba baca dari environment variable dulu (untuk production/Railway)
firebase_creds_json = os.getenv("FIREBASE_CREDENTIALS")

if firebase_creds_json:
    # Di server Railway: kredensial dari environment variable
    cred_dict = json.loads(firebase_creds_json)
    cred = credentials.Certificate(cred_dict)
else:
    # Di lokal: tetap baca dari file seperti biasa
    cred_path = os.path.join(os.path.dirname(__file__), "..", "firebase-service-account.json")
    cred = credentials.Certificate(cred_path)

firebase_app = firebase_admin.initialize_app(cred)
db = firestore.client()