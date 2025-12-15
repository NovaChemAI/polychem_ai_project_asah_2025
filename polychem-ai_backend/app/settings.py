import os

# =========================================================
# Koyeb-safe paths
# /app  -> READ ONLY
# /tmp  -> WRITABLE
# =========================================================

# Root static directory (writable)
STATIC_DIR = os.getenv("STATIC_DIR", "/tmp/static")

# Folder simpan gambar compounds
COMPOUNDS_DIR = os.path.join(STATIC_DIR, "compounds")