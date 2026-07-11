from fastapi import Header, HTTPException
from firebase_admin import auth
from .firebase_admin_config import firebase_app
from typing import Optional
from fastapi import Header


async def get_current_user(authorization: str = Header(...)):
    if not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Token tidak valid")
    
    id_token = authorization.split("Bearer ")[1]
    
    try:
        decoded_token = auth.verify_id_token(id_token)
        return decoded_token
    except Exception:
        raise HTTPException(status_code=401, detail="Token tidak valid atau kedaluwarsa")
    
    
async def get_optional_user(authorization: Optional[str] = Header(None)):
    if not authorization or not authorization.startswith("Bearer "):
        return None

    id_token = authorization.split("Bearer ")[1]
    try:
        return auth.verify_id_token(id_token)
    except Exception:
        return None