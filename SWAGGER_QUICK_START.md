# ⚡ Swagger Testing - Super Quick Guide

**Goal:** Test backend manually in browser (no coding needed!)

---

## 🎯 3 Steps Only

### Step 1: Start Backend (Terminal 1)
```bash
cd polychem-ai_backend
python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```
✅ Wait for: `"Dataset loaded!"`

### Step 2: Open Browser
```
URL: http://127.0.0.1:8000/docs
```

### Step 3: Test Endpoints
Click endpoint → "Try it out" → Enter data → "Execute" → See results!

---

## 🧪 Quick Test (2 minutes)

### Test 1: Health Check
```
1. Find: GET /health (green)
2. Click: "Try it out"
3. Click: "Execute"
4. See: {"status": "ok"} ✅
```

### Test 2: Predict (Main Test)
```
1. Find: POST /predict (blue)
2. Click: "Try it out"
3. Enter:
   {
     "smiles": "CCO"
   }
4. Click: "Execute"
5. See: Predictions + similar compounds ✅
```

### Test 3: History
```
1. Find: GET /history (green)
2. Click: "Try it out"
3. Click: "Execute"
4. See: Your past tests ✅
```

---

## 🧬 SMILES to Test

```
Simple:
  "C"              - Methane
  "CC"             - Ethane
  "CCO"            - Ethanol ← Best first test

Complex:
  "c1ccccc1"       - Benzene
  "CC(=O)O"        - Acetic acid
  "CC(=O)Oc1ccccc1C(=O)O"  - Aspirin

Invalid (tests error handling):
  "invalid"        - Should get error
  "###"            - Should get error
```

---

## 📊 What You'll See

**Response includes:**
```
✅ Predicted compound name
✅ Chemical formula (e.g., C2H6O)
✅ Glass Transition Temp (Tg)
✅ Image URL (2D structure)
✅ 3 Similar compounds
✅ Why they're similar
```

---

## ⏱️ Speed

- **First request:** 5-20 seconds ⏳
- **Second request (cached):** 50 ms ⚡
- **Health check:** <1 sec ✨

---

## 🎯 Test Different SMILES

| SMILES | Description | Try It |
|--------|-------------|--------|
| `"C"` | Simplest | First |
| `"CC"` | 2 carbons | Second |
| `"CCO"` | Ethanol | Third ← Recommended |
| `"c1ccccc1"` | Benzene ring | Advanced |

---

## 💡 Pro Tips

1. **Test simple first** - C, CC, CCO
2. **Watch the time** - Notice fast caching
3. **Try invalid** - See error handling (e.g., "invalid")
4. **Check history** - See all your requests

---

## 🔗 URLs

| Purpose | URL |
|---------|-----|
| **Swagger UI** | http://127.0.0.1:8000/docs |
| **Alternative Docs** | http://127.0.0.1:8000/redoc |
| **Health Check** | http://127.0.0.1:8000/health |

---

## ❌ If Backend Doesn't Start

```
Error: "Cannot GET /docs"
Fix: Run backend first!
  python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```

---

## ✅ Success = This Response

```json
{
  "status": "success",
  "input_smiles": "CCO",
  "new_compound": {
    "name": "Ethanol-derived Polyether",
    "formula": "C2H6O",
    "molecular_weight": 46.04,
    "tg": 45.5,
    "justifikasi": "Simple alcohol-ether structure..."
  },
  "similar_compounds": [
    {
      "rank": 1,
      "similarity_percent": 95.2,
      "justifikasi": "Same functional group..."
    },
    ...
  ]
}
```

---

**Status:** ✅ Ready to test!

**Go to:** http://127.0.0.1:8000/docs
