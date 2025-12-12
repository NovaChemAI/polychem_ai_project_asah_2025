import os

# BASE_DIR = root project (novachem-ai-backend)
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

# folder static 
STATIC_DIR = os.path.join(BASE_DIR, "static")

# folder simpan gambar compunds
COMPOUND_DIR = os.path.join(STATIC_DIR, "compounds")
