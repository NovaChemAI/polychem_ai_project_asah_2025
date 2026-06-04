# 🔧 Backend - Critical Fixes Implementation Guide

**Status:** 3 Critical Fixes Needed  
**Estimated Time:** 30 minutes  
**Difficulty:** Easy (mostly 1-2 line changes)

---

## Fix #1: CORS Too Permissive 🔴 CRITICAL

### ❌ Current Code
**File:** `app/main.py` (lines 47-63)

```python
# Current code accepts requests from ANYWHERE
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # ⚠️ SECURITY RISK!
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

**Why This Is Bad:**
- Allows CSRF attacks from any domain
- Exposes API to misuse
- Not secure for production

### ✅ Fixed Code

```python
# app/main.py - Replace the CORS section (lines 47-63)

# Get frontend URL from environment (production) or default to localhost (dev)
frontend_url = os.getenv("FRONTEND_URL", "http://localhost:5173")

# Whitelist only specific origins
allow_origins = [
    frontend_url,              # Production frontend
    "http://localhost:5173",   # Local Vite dev
    "http://127.0.0.1:5173",   # Alternative localhost
    "http://localhost:3000",   # Alternative port
    "http://127.0.0.1:3000",   # Alternative localhost
]

# Remove duplicates
allow_origins = list(set(allow_origins))

app.add_middleware(
    CORSMiddleware,
    allow_origins=allow_origins,  # ✅ Specific whitelist
    allow_credentials=False,
    allow_methods=["GET", "POST"],  # Only needed methods
    allow_headers=["Content-Type"],  # Only needed headers
)
```

### Complete Updated File
**Replace lines 47-63 in `app/main.py` with:**

```python
# =========================================================
# CORS Configuration (FIXED: Specific whitelist)
# =========================================================
# Get frontend URL from environment (production) or default to localhost (dev)
frontend_url = os.getenv("FRONTEND_URL", "http://localhost:5173")

# Build whitelist of allowed origins
allow_origins = [
    frontend_url,              # Production frontend
    "http://localhost:5173",   # Local Vite dev
    "http://127.0.0.1:5173",   # Alternative localhost
    "http://localhost:3000",   # Alternative port
    "http://127.0.0.1:3000",   # Alternative localhost
]

# Remove duplicates in case frontend_url matches a dev URL
allow_origins = list(set(allow_origins))

app.add_middleware(
    CORSMiddleware,
    allow_origins=allow_origins,  # ✅ FIXED: Specific whitelist instead of ["*"]
    allow_credentials=False,
    allow_methods=["GET", "POST"],  # ✅ FIXED: Only needed methods
    allow_headers=["Content-Type"],  # ✅ FIXED: Only needed headers
)
```

### For Koyeb Production
**Add to Koyeb Environment Variables:**
```
FRONTEND_URL=https://your-frontend-domain.vercel.app
```

---

## Fix #2: API Key Exposed 🔴 CRITICAL

### ❌ Current Problem
**File:** `.env` (in repository)

```env
GOOGLE_API_KEY=AIzaSyDgcPbMWYYkyj4pkyQX474LcuDhU5w_iNU  # ⚠️ EXPOSED!
```

**Why This Is Bad:**
- API key visible in git history forever
- Anyone with repo access can abuse API
- API key not rotated
- Anyone finding this key can make unlimited API calls

### ✅ Solution: Rotate Key

**Step 1: Create New API Key**
1. Go to [Google Cloud Console](https://console.cloud.google.com)
2. Select your project
3. Go to "APIs & Services" → "Credentials"
4. Click "Create Credentials" → "API Key"
5. Copy the new key

**Step 2: Update .env**
```env
GOOGLE_API_KEY=your_new_api_key_here
ID_DATASET_DRIVE=13L5ZFx_vZyrTwS4tUeWADO1-SE_oPHkM
```

**Step 3: Disable Old Key**
1. In Google Cloud Console
2. Find the old key
3. Click to select it
4. Click "Delete" or "Disable"

**Step 4: Create .env.example Template**

Create file `polychem-ai_backend/.env.example`:
```env
# Google Generative AI
GOOGLE_API_KEY=your_api_key_here

# Google Drive Dataset ID
ID_DATASET_DRIVE=13L5ZFx_vZyrTwS4tUeWADO1-SE_oPHkM

# Optional Configuration
FRONTEND_URL=http://localhost:5173
STATIC_DIR=/tmp/static
DATA_CACHE_DIR=/tmp/data_cache
CACHE_DIR=/tmp/polychem_cache
HISTORY_LIMIT=10
HISTORY_TTL_SECONDS=3600
```

**Step 5: Add .env to .gitignore**
```bash
# Run this in terminal
echo ".env" >> .gitignore
echo ".env.local" >> .gitignore
```

**Step 6: Update Koyeb**
1. Go to Koyeb dashboard
2. Select your deployment
3. Go to Settings → Environment
4. Update `GOOGLE_API_KEY` with new key
5. Deploy

---

## Fix #3: Input Validation Missing 🔴 CRITICAL

### ❌ Current Code
**File:** `app/schemas.py` (line 4)

```python
class RecommendRequest(BaseModel):
    smiles: str = Field(..., min_length=1)  # ⚠️ No max_length!
```

**Why This Is Bad:**
- User could send 1 MB string → memory DoS
- RDKit would try to parse huge string → slow
- Server resource exhaustion possible

### ✅ Fixed Code
**Replace line 4 in `app/schemas.py` with:**

```python
class RecommendRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)  # ✅ FIXED: Added max_length
```

**Why 500?**
- Realistic SMILES rarely exceed 200 characters
- 500 is safe buffer for complex molecules
- Prevents abuse/DoS attempts

---

## 📝 Summary of Changes

| File | Line | Change | Time |
|------|------|--------|------|
| `app/main.py` | 47-63 | CORS: `["*"]` → specific whitelist | 2 min |
| `.env` | All | Rotate API key | 5 min |
| `.env.example` | New | Create template | 2 min |
| `.gitignore` | Add | Add `.env` | 1 min |
| `app/schemas.py` | 4 | Add `max_length=500` | 1 min |
| `Koyeb dashboard` | N/A | Update FRONTEND_URL env var | 5 min |
| **TOTAL** | | **6 changes** | **~30 min** |

---

## Step-by-Step Implementation

### Step 1: Fix CORS in app/main.py
1. Open `polychem-ai_backend/app/main.py`
2. Find the CORS section (around line 47)
3. Replace `allow_origins=["*"]` with the fixed code
4. Save file

### Step 2: Rotate API Key
1. Go to Google Cloud Console
2. Create new API key
3. Update `.env` with new key
4. Mark old key as disabled in Google Console
5. Save file

### Step 3: Add .env.example
1. Create new file `polychem-ai_backend/.env.example`
2. Copy the template content
3. Replace actual values with placeholders
4. Save file

### Step 4: Update .gitignore
1. Open `polychem-ai_backend/.gitignore` (create if doesn't exist)
2. Add:
   ```
   .env
   .env.local
   ```
3. Save file

### Step 5: Fix Input Validation
1. Open `polychem-ai_backend/app/schemas.py`
2. Find `RecommendRequest` class
3. Change `min_length=1` to `min_length=1, max_length=500`
4. Save file

### Step 6: Update Koyeb
1. Go to [Koyeb Dashboard](https://app.koyeb.com)
2. Select your deployment
3. Click "Settings" → "Environment"
4. Update:
   - `GOOGLE_API_KEY=` (new key)
   - `FRONTEND_URL=https://your-frontend.vercel.app`
5. Click "Update" and redeploy

---

## 🧪 Verification

After making changes:

```bash
# 1. Test locally
cd polychem-ai_backend
pip install -r requirements.txt
python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000

# 2. In another terminal, run tests
python test_backend.py

# 3. Should see "8/8 tests passed"
```

### Expected Test Output
```
============================================================
                 PolyChem AI Backend Test Suite
============================================================
Target URL: http://127.0.0.1:8000

✓ PASS | Health Check
✓ PASS | Predict: Methane
✓ PASS | Predict: Ethane
✓ PASS | Predict: Ethanol
✓ PASS | Predict: Benzene
✓ PASS | Handle invalid input
✓ PASS | Get history
✓ PASS | Concurrent requests

============================================================
Total: 8/8 tests passed
All tests passed!
============================================================
```

---

## ✅ Verification Checklist

- [ ] Changed CORS `allow_origins=["*"]` to specific whitelist
- [ ] Created new Google API key
- [ ] Updated `.env` with new key
- [ ] Marked old API key as disabled in Google Console
- [ ] Created `.env.example` template
- [ ] Added `.env` to `.gitignore`
- [ ] Added `max_length=500` to SMILES input validation
- [ ] Tested locally (8/8 tests pass)
- [ ] Updated Koyeb environment variables
- [ ] Redeployed to Koyeb
- [ ] Health check on production: `curl https://slim-danika-polychem-ab276767.koyeb.app/health`

---

## ⚠️ Important Notes

1. **Old API Key:** Once rotated, the old key should never be used again. Disable it in Google Cloud.

2. **Git History:** The exposed key is still in git history. Consider:
   - Using git filter to remove from history, OR
   - Creating new repo if this is early stage

3. **FRONTEND_URL:** 
   - Local dev: `http://localhost:5173`
   - Production: `https://your-actual-frontend.vercel.app`

4. **Testing:** After changes, run full test suite to ensure nothing broke.

---

## 📞 Questions?

If any step is unclear or tests fail, the troubleshooting guide is in `BACKEND_SETUP_READY.md`.

**Status After Fixes:** Backend will be 95% production-ready ✅
