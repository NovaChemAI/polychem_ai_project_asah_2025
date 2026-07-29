from fastapi import Header, HTTPException
from firebase_admin import auth
from .firebase_admin_config import firebase_app
from typing import Optional
from fastapi import Header


async def get_current_user(authorization: str = Header(...)):
    """
    Dipakai untuk endpoint yang WAJIB email sudah terverifikasi
    (misal: /predict, /library, /history/mine).
    """
    if not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Token tidak valid")

    id_token = authorization.split("Bearer ")[1]

    try:
        decoded_token = auth.verify_id_token(id_token)
    except Exception:
        raise HTTPException(status_code=401, detail="Token tidak valid atau kedaluwarsa")

    if not decoded_token.get("email_verified", False):
        raise HTTPException(status_code=403, detail="Email belum diverifikasi. Silakan cek inbox kamu.")

    return decoded_token


async def get_current_user_no_verify(authorization: str = Header(...)):
    """
    Sama seperti get_current_user, tapi TIDAK mengecek email_verified.
    Dipakai khusus untuk endpoint yang perlu diakses SEBELUM
    user verifikasi email (misal: /auth/sync-profile saat baru daftar).
    """
    if not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Token tidak valid")

    id_token = authorization.split("Bearer ")[1]

    try:
        decoded_token = auth.verify_id_token(id_token)
    except Exception:
        raise HTTPException(status_code=401, detail="Token tidak valid atau kedaluwarsa")

    return decoded_token


async def get_optional_user(authorization: Optional[str] = Header(None)):
    if not authorization or not authorization.startswith("Bearer "):
        return None

    id_token = authorization.split("Bearer ")[1]
    try:
        return auth.verify_id_token(id_token)
    except Exception:
        return None