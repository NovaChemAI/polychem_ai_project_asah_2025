# 🧪 Backend SMILES Testing Guide

**Purpose:** Test the backend API with real SMILES chemical notation  
**Date:** January 28, 2026

---

## 📋 What is SMILES?

SMILES = Simplified Molecular Input Line Entry System

**Examples:**
```
C           = Methane (CH4)
CC          = Ethane (C2H6)
CCO         = Ethanol (C2H5OH)
c1ccccc1    = Benzene (C6H6)
CC(=O)O     = Acetic Acid (CH3COOH)
```

**Why?** Chemical notation easy to input and process

---

## 🚀 How to Test

### Step 1: Start Backend Server

**Terminal 1:**
```bash
cd polychem-ai_backend
pip install -r requirements.txt
python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```

Wait for: `"Dataset loaded!"` message

### Step 2: Test with Script

**Terminal 2:**
```bash
cd polychem-ai_backend
python test_smiles_interactive.py
```

### Step 3: View Results

Script will test:
- ✅ 8 valid SMILES
- ✅ 3 invalid SMILES (should be rejected)
- ✅ Response times
- ✅ Prediction accuracy

---

## 📊 Test SMILES Strings

### Valid SMILES (Should Work)

| SMILES | Name | Description |
|--------|------|-------------|
| `C` | Methane | Simplest organic compound |
| `CC` | Ethane | 2 carbon atoms |
| `CCO` | Ethanol | Common alcohol |
| `CCOC` | Ethyl Methyl Ether | Ether compound |
| `c1ccccc1` | Benzene | Aromatic ring |
| `CC(=O)O` | Acetic Acid | Carboxylic acid |
| `CC(=O)Oc1ccccc1C(=O)O` | Aspirin | Well-known drug |
| `CCCCCCCCCCCCc1ccccc1` | Dodecylbenzene | Long chain |

### Invalid SMILES (Should Reject)

| SMILES | Reason |
|--------|--------|
| `random string` | Not chemical notation |
| `###` | Invalid characters |
| `C@@@C` | Invalid bond notation |

---

## 💡 What the Script Tests

### 1. Health Check
```
✓ Verifies backend is running
✓ Expected: {"status": "ok"}
```

### 2. Valid SMILES Tests
```
For each SMILES:
  ✓ Send to /predict endpoint
  ✓ Parse response
  ✓ Show predictions
  ✓ Show similar compounds
  ✓ Measure response time
```

### 3. Invalid SMILES Tests
```
For invalid inputs:
  ✓ API should reject (400 error)
  ✓ Show error message
  ✓ Verify graceful error handling
```

---

## 📈 Expected Output

```
================================================================================
              POLYCHEM AI - SMILES PREDICTION TESTING
================================================================================

Checking backend health...

✓ Backend is running on http://127.0.0.1:8000
✓ Health check: {'status': 'ok'}

================================================================================
                        VALID SMILES TESTS
================================================================================

────────────────────────────────────────────────────────────
Testing: CCO
Description: Ethanol - common alcohol
────────────────────────────────────────────────────────────
Response Time: 2.34s
Status Code: 200
✓ Prediction successful!

New Compound Prediction:
  Name: Ethanol-derived Polyether
  Formula: C2H6O
  Molecular Weight: 46.04 g/mol
  Tg (Glass Transition): 45.5°C
  Justification: Simple ether linkage with unique alcohol functionality
  Image: /static/compounds/compound_xyz.png

Similar Compounds (Top 3):

  #1 Similarity: 95.2%
     SMILES: CCOC
     Name: Ethyl methyl ether
     Tg: 42.1°C
     Reason: Same ether functional group

  #2 Similarity: 89.3%
     SMILES: CC(C)O
     Name: Isopropanol
     Tg: 38.7°C
     Reason: Similar alcohol structure

  #3 Similarity: 82.1%
     SMILES: CCN
     Name: Ethylamine
     Tg: 35.2°C
     Reason: Similar chain length and atom types

... (more tests)

================================================================================
                    INVALID SMILES TESTS
================================================================================

────────────────────────────────────────────────────────────
Testing Invalid SMILES: ###
────────────────────────────────────────────────────────────
✓ Correctly rejected invalid SMILES
Status: 400
Error: Invalid SMILES: ### (invalid characters)

... (more tests)

================================================================================
                        TEST SUMMARY
================================================================================
Valid SMILES: 8/8 passed
Invalid SMILES: 3/3 passed

Total: 11/11 tests passed
✓ All tests passed!
```

---

## 🔍 Response Details Explained

### New Compound Fields
```json
{
  "name": "Ethanol-derived Polyether",
  "formula": "C2H6O",
  "molecular_weight": 46.04,
  "tg": 45.5,                    // Glass Transition Temperature in °C
  "justifikasi": "Simple ether...",
  "image_url": "/static/compounds/..."
}
```

### Similar Compounds Fields
```json
{
  "rank": 1,                      // 1st, 2nd, or 3rd most similar
  "smiles": "CCOC",
  "similarity_score": 0.952,      // 0-1 scale
  "similarity_percent": 95.2,     // 0-100 scale
  "justifikasi": "Same ether..."
}
```

---

## ⚡ Performance Expectations

### First Request (New SMILES)
- Dataset load: 1-30 sec (one-time)
- Processing: 5-20 sec (depends on Google API)
- Total: 5-20 sec

### Subsequent Requests (Cached)
- L1 Cache hit: ~50 ms (instant!)
- Total: <100 ms

### Response Times by Phase
```
Fingerprint generation:  0.1 sec
Similarity search:       0.05 sec
LLM calls:              5-15 sec  ← Slowest (Google Gemini)
Image generation:       0.2 sec
Total:                  5-20 sec
```

---

## 🐛 Troubleshooting

### Error: "Cannot connect to backend"
```
Solution: Make sure backend is running
$ python -m uvicorn app.main:app --reload --host 127.0.0.1 --port 8000
```

### Error: "Request timeout"
```
Solution: First run takes time for dataset download
- Wait 30+ seconds on first request
- Subsequent requests will be faster
```

### Error: "Invalid SMILES"
```
Solution: Check SMILES notation
- Valid: C, CC, CCO, c1ccccc1
- Invalid: random text, ###, etc.
```

### Error: "Module not found"
```
Solution: Install dependencies
$ pip install -r requirements.txt
```

---

## 🎯 Manual Testing (curl)

If you prefer testing with curl instead:

### Test Health
```bash
curl http://127.0.0.1:8000/health
```

### Test Prediction
```bash
curl -X POST http://127.0.0.1:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

### Get History
```bash
curl http://127.0.0.1:8000/history
```

---

## 📊 Test Results Checklist

After running tests, verify:

- [ ] Health check passes
- [ ] All 8 valid SMILES tests pass
- [ ] All 3 invalid SMILES tests rejected
- [ ] Response times < 30 seconds
- [ ] Predictions make chemical sense
- [ ] Similar compounds are actually similar
- [ ] Images are generated successfully
- [ ] No errors in server terminal

---

## 🎓 Learning the Results

### How to Read Predictions

**Tg (Glass Transition Temperature):**
- Measured in °C
- Range typically: -80°C to 250°C
- Lower Tg = more flexible polymer
- Higher Tg = more rigid polymer

**Similarity Score:**
- 0.0 = completely different
- 0.5 = somewhat similar
- 1.0 = identical
- Calculated using Tanimoto coefficient

**Justification:**
- Why is this compound interesting?
- Why are similar compounds similar?
- Generated by AI (with fallback heuristics)

---

## 📝 Common SMILES for Testing

**Small Molecules:**
```
C             Methane
CC            Ethane
CCC           Propane
CCCC          Butane
CCO           Ethanol
CC(C)C        Isobutane
```

**Aromatic:**
```
c1ccccc1      Benzene
Cc1ccccc1     Toluene
c1ccccc1c1ccccc1  Biphenyl
```

**Functional Groups:**
```
CC(=O)O       Carboxylic acid (acetic acid)
CCN           Amine
CCOC          Ether
CC(=O)N       Amide
```

---

## 🚀 Next Steps

1. **Run the test script**
   ```bash
   python test_smiles_interactive.py
   ```

2. **Review results**
   - Check if all tests pass
   - Note response times
   - Verify predictions make sense

3. **If issues found**
   - Check backend logs
   - Read troubleshooting section
   - Refer to BACKEND_COMPREHENSIVE_REPORT.md

4. **Manual testing**
   - Try different SMILES strings
   - Test edge cases
   - Verify error handling

---

## 📞 Questions?

**Q: Why is first request slow?**  
A: Must download 7284 SMILES from Google Drive. Cached after.

**Q: Can I test with any SMILES?**  
A: Yes! As long as it's valid chemical notation.

**Q: How many requests can I make?**  
A: Limited by backend rate limit (10/min per IP).

**Q: Can I test with file upload?**  
A: Not yet - API only accepts SMILES string parameter.

---

**Status:** Ready to test!

**Next Step:** Run `python test_smiles_interactive.py` and see results
