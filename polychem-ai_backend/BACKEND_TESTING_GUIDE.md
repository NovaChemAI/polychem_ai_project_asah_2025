# Backend Testing Guide - PolyChem AI

Panduan lengkap untuk setup, menjalankan, dan test Backend FastAPI.

## Prerequisites

Pastikan Anda sudah install:
- Python 3.8+ 
- pip atau conda
- Git

## Quick Start

### 1. Install Dependencies

```bash
cd polychem-ai_backend
pip install -r requirements.txt
```

Jika ada error dengan `rdkit`, coba:
```bash
pip install rdkit-pypi
# atau untuk conda users:
conda install -c conda-forge rdkit
```

### 2. Setup Environment Variables

File `.env` sudah ada dengan:
```
GOOGLE_API_KEY=AIzaSyDgcPbMWYYkyj4pkyQX474LcuDhU5w_iNU
```

Pastikan API key valid, atau tambahkan di `.env`:
```
GOOGLE_API_KEY=your_api_key_here
FRONTEND_URL=http://localhost:5173  # untuk development
```

### 3. Run Backend Server

**Development Mode (with auto-reload):**
```bash
python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```

**Production Mode:**
```bash
python -m uvicorn app.main:app --host 0.0.0.0 --port 8000
```

Anda akan melihat output:
```
INFO:     Uvicorn running on http://127.0.0.1:8000
INFO:     Application startup complete
Loading dataset...
Dataset loaded!
```

### 4. Verify Backend is Running

Buka di browser atau terminal:
```bash
curl http://127.0.0.1:8000/health
# Expected response: {"status":"ok"}
```

---

## Testing

### Option A: Menggunakan Test Script (Recommended)

Kami sudah siapkan comprehensive test suite: `test_backend.py`

```bash
# Terminal 1: Start Backend
python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000

# Terminal 2: Run Tests
python test_backend.py
```

**Test Suite includes:**
- ✓ Health check
- ✓ Valid SMILES prediction
- ✓ Invalid SMILES handling
- ✓ History endpoint
- ✓ Response format validation
- ✓ Performance testing
- ✓ Sequential requests
- ✓ Edge cases

### Option B: Manual Testing

#### 1. Health Check
```bash
curl http://127.0.0.1:8000/health
```

#### 2. Predict with SMILES
```bash
curl -X POST http://127.0.0.1:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

#### 3. Get History
```bash
curl http://127.0.0.1:8000/history
```

#### 4. Interactive Testing (Python)
```python
import requests

# Simple predict
response = requests.post(
    "http://127.0.0.1:8000/predict",
    json={"smiles": "CCO"}
)
print(response.json())
```

### Option C: Using Postman

Import ini ke Postman:

**Health Check**
- Method: GET
- URL: `http://127.0.0.1:8000/health`

**Predict**
- Method: POST
- URL: `http://127.0.0.1:8000/predict`
- Body (JSON):
  ```json
  {
    "smiles": "CCO"
  }
  ```

**History**
- Method: GET
- URL: `http://127.0.0.1:8000/history`

---

## API Endpoints

### 1. POST /predict
**Request:**
```json
{
  "smiles": "CCO"
}
```

**Response:**
```json
{
  "status": "success",
  "input_smiles": "CCO",
  "new_compound": {
    "name": "...",
    "smiles": "...",
    "formula": "...",
    "molecular_weight": 46.04,
    "tg": 100.5,
    "tg_justification": "...",
    "polymer_class": "...",
    "justifikasi": "...",
    "image_url": "/static/compounds/compound_abc123.png",
    "fingerprint_length": 2048,
    "image_filename": "compound_abc123.png"
  },
  "similar_compounds": [
    {
      "rank": 1,
      "smiles": "...",
      "name": "...",
      "formula": "...",
      "molecular_weight": 0.0,
      "tg": 0.0,
      "similarity_score": 0.95,
      "similarity_percent": 95.0,
      "justifikasi": "...",
      "image_url": "...",
      "image_filename": "..."
    }
  ]
}
```

### 2. GET /health
**Response:** `{"status": "ok"}`

### 3. GET /history
**Response:** Array of past predictions

---

## Common Issues & Solutions

### Issue 1: "Module not found: rdkit"
**Solution:**
```bash
# Uninstall and reinstall
pip uninstall rdkit
pip install rdkit-pypi

# Or via conda
conda install -c conda-forge rdkit
```

### Issue 2: "GOOGLE_API_KEY not set"
**Solution:**
- Pastikan file `.env` ada di folder backend
- Check bahwa API key valid di Google Cloud Console
- Restart server setelah mengubah `.env`

### Issue 3: "Port 8000 already in use"
**Solution:**
```bash
# Cari proses yang pakai port 8000
lsof -i :8000  # Linux/Mac
netstat -ano | findstr :8000  # Windows

# Kill process atau gunakan port berbeda
python -m uvicorn app.main:app --reload --port 8001
```

### Issue 4: "Dataset not loading" atau timeout pada /predict
**Solution:**
- Download dataset dari Google Drive mungkin lambat di first run
- Cek connection, atau download manual dan edit `app/store.py`
- Increase timeout: `timeout=60` di test script

### Issue 5: "Image generation failed"
**Solution:**
- Pastikan `/tmp/static/compounds` writable
- Cek permission folder
- Atau set `STATIC_DIR` environment variable

---

## Performance Tuning

### 1. Dataset Caching
Dataset di-cache di `/tmp/data_cache/dataset.csv` setelah download pertama.
Hapus untuk force re-download:
```bash
rm -rf /tmp/data_cache
```

### 2. Fingerprint Caching
Cache L2 diskcache di `/tmp/polychem_cache`. Hapus untuk reset:
```bash
rm -rf /tmp/polychem_cache
```

### 3. History Cleanup
History auto-cleanup berdasarkan TTL (default 1 jam). Ubah via env:
```bash
export HISTORY_TTL_SECONDS=3600  # 1 jam
export HISTORY_LIMIT=20          # Simpan max 20 items
```

---

## Environment Variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `GOOGLE_API_KEY` | - | Google Generative AI API Key (REQUIRED) |
| `FRONTEND_URL` | (auto) | Frontend URL untuk CORS |
| `STATIC_DIR` | `/tmp/static` | Path untuk simpan gambar |
| `DATA_CACHE_DIR` | `/tmp/data_cache` | Path untuk cache dataset |
| `CACHE_DIR` | `/tmp/polychem_cache` | Path untuk diskcache |
| `CACHE_VERSION` | `v1` | Cache version (ubah untuk invalidate) |
| `HISTORY_LIMIT` | `10` | Max history items |
| `HISTORY_TTL_SECONDS` | `3600` | History time-to-live |
| `ID_DATASET_DRIVE` | (set) | Google Drive file ID untuk dataset |
| `GEMINI_MODEL_FAST` | `gemini-2.5-flash` | Model untuk naming |
| `GEMINI_MODEL_TG` | `gemini-2.5-flash` | Model untuk Tg prediction |

---

## Troubleshooting Steps

1. **Check Server is Running**
   ```bash
   curl http://127.0.0.1:8000/health
   ```

2. **Check Logs**
   - Lihat output di terminal tempat uvicorn running
   - Error biasanya di-print langsung

3. **Test Simple Input**
   ```bash
   curl -X POST http://127.0.0.1:8000/predict \
     -H "Content-Type: application/json" \
     -d '{"smiles": "C"}'
   ```

4. **Check Network**
   ```bash
   # Windows
   netstat -ano | findstr :8000
   
   # Linux/Mac
   lsof -i :8000
   ```

5. **Restart Backend**
   - Kill process
   - Clear cache: `rm -rf /tmp/polychem_*`
   - Restart: `python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000`

---

## Next Steps

1. ✓ Backend running dan tested
2. [ ] Test Frontend integration
3. [ ] Check CORS configuration
4. [ ] Test Firebase authentication
5. [ ] Deploy ke production (Koyeb)

---

Generated: January 28, 2026
