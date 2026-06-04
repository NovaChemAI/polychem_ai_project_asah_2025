# ✅ Manual Testing Checklist - Polychem AI Project

**Test Date:** ******\_\_\_******  
**Tester:** Dikky  
**Environment:** Development (localhost)

---

## 🔧 Pre-Flight Checks

### Backend Status

- [ ] Backend running on `http://localhost:8000`
- [ ] `/health` endpoint returns `{"status": "healthy"}`
- [ ] Swagger UI accessible at `http://localhost:8000/docs`

### Frontend Status

- [ ] Frontend running on `http://localhost:5173`
- [ ] No console errors related to Firebase config
- [ ] Hard refresh (Ctrl+F5) loads properly without white screen

---

## 🧪 Critical Path Testing

### 1. Firebase Authentication

**Priority:** HIGH 🔴  
**Location:** LoginPage, RegisterPage

#### Test 1.1: User Registration

1. Navigate to `http://localhost:5173/register`
2. Fill form:
   - **Name:** Test User
   - **Email:** testuser@example.com
   - **Password:** TestPass123!
3. Click "Register"
4. **Expected:** Toast success ✅ + redirect to `/`

#### Test 1.2: User Login

1. Navigate to `http://localhost:5173/login`
2. Fill credentials from Test 1.1
3. Click "Login"
4. **Expected:** Toast success ✅ + redirect to `/`

#### Test 1.3: Logout

1. Click user profile/logout button in sidebar
2. **Expected:** Redirect to `/login`

#### Test 1.4: Firebase Console Validation

1. Open Firebase Console → Authentication
2. **Expected:** See `testuser@example.com` in user list

---

### 2. API Integration (Backend ↔ Frontend)

**Priority:** HIGH 🔴  
**Location:** HomePage `/predict` endpoint

#### Test 2.1: Valid SMILES Prediction

1. Navigate to `http://localhost:5173/` (logged in)
2. Input SMILES: `CC(C)(C)c1ccc(O)cc1` (BPA monomer)
3. Click "Generate Prediction"
4. **Expected:**
   - Loading spinner shows
   - Redirects to `/chemical-detail/:id`
   - Displays:
     - ✅ New compound name
     - ✅ Molecular structure PNG image
     - ✅ Properties table (molecular_weight, LogP, TPSA)
     - ✅ AI description paragraph
     - ✅ Similar compounds section (3-5 items)

#### Test 2.2: Invalid SMILES Error Handling

1. Navigate to `http://localhost:5173/`
2. Input invalid: `XXXXXX`
3. Click "Generate Prediction"
4. **Expected:**
   - Error toast shows ❌
   - Stays on HomePage (no redirect)
   - Backend returns 400 status

#### Test 2.3: Empty Input Validation

1. Navigate to `http://localhost:5173/`
2. Leave input empty, click "Generate Prediction"
3. **Expected:** Frontend validation error (no API call)

---

### 3. Data Persistence (Firestore)

**Priority:** MEDIUM 🟡  
**Location:** LibraryPage, HistoryPage

#### Test 3.1: Save to Library

1. After Test 2.1 (on ChemicalDetailPage)
2. Click "Save to Library" button
3. **Expected:**
   - Toast success "Saved to library" ✅
   - Firebase Firestore → `predictions` collection gets new doc

#### Test 3.2: View Library

1. Navigate to `http://localhost:5173/library`
2. **Expected:**
   - Shows saved prediction from Test 3.1
   - Card displays compound name + properties
   - Click card → navigates to detail page

#### Test 3.3: View History

1. Navigate to `http://localhost:5173/history`
2. **Expected:**
   - Shows all predictions made (Test 2.1)
   - Chronological order (newest first)
   - Each entry clickable → detail page

---

### 4. Static Assets & Image Serving

**Priority:** MEDIUM 🟡  
**Location:** ChemicalDetailPage

#### Test 4.1: Molecular Structure Images

1. After Test 2.1, inspect ChemicalDetailPage
2. Open DevTools → Network tab, filter Images
3. **Expected:**
   - PNG image loads from `http://localhost:8000/static/compounds/[hash].png`
   - Image visible (not broken icon)
   - Similar compounds also show images

#### Test 4.2: Image Caching

1. Repeat Test 2.1 with same SMILES (BPA)
2. Check backend terminal logs
3. **Expected:**
   - Image served from cache (no regeneration log)
   - Faster response time (~100ms vs initial generation)

---

### 5. Error Recovery & Edge Cases

**Priority:** MEDIUM 🟡

#### Test 5.1: Backend Offline Recovery

1. Stop backend terminal (Ctrl+C)
2. Try Test 2.1 (generate prediction)
3. **Expected:**
   - Frontend shows error toast: "Failed to fetch" or network error ❌
   - No crash, stays on HomePage
4. Restart backend, retry
5. **Expected:** Works normally ✅

#### Test 5.2: Missing Image Fallback

1. Manually delete an image from `polychem-ai_backend/static/compounds/`
2. Navigate to that prediction's detail page
3. **Expected:**
   - Shows placeholder or broken image gracefully
   - No console errors crash the page

#### Test 5.3: Large Dataset Handling

1. Input complex polymer SMILES: `C(C(C(C(C(C(C(C(Cl)Cl)Cl)Cl)Cl)Cl)Cl)Cl` (Octachlorobutane)
2. **Expected:**
   - Backend may take 3-5s (LLM generation)
   - No timeout errors
   - Returns valid result

---

### 6. UI/UX & Responsiveness

**Priority:** LOW 🟢

#### Test 6.1: Dark Mode Toggle

1. Click theme toggle in navbar/sidebar
2. **Expected:**
   - Entire UI switches to dark mode
   - LocalStorage saves preference
   - Refresh retains dark mode

#### Test 6.2: Mobile View (DevTools)

1. Open DevTools → Device Toolbar (Ctrl+Shift+M)
2. Select "iPhone 12 Pro"
3. Navigate through all pages
4. **Expected:**
   - Sidebar collapses to hamburger menu
   - Cards stack vertically
   - Text remains readable

#### Test 6.3: Toast Notifications

1. Trigger success (Test 3.1) and error (Test 2.2) toasts
2. **Expected:**
   - Appears in top-right/center
   - Auto-dismiss after 3-5s
   - Distinct styling (green vs red)

---

## 🔬 Advanced Testing (Optional)

### 7. Backend API Direct Testing

**Tool:** Postman, cURL, or Swagger UI

#### Test 7.1: `/predict` POST

```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(C)(C)c1ccc(O)cc1"}'
```

**Expected:** JSON with `new_compound`, `similar_compounds`, `prediction_id`

#### Test 7.2: `/health` GET

```bash
curl http://localhost:8000/health
```

**Expected:** `{"status": "healthy"}`

#### Test 7.3: `/history` GET

```bash
curl http://localhost:8000/history
```

**Expected:** Array of past predictions

---

### 8. Machine Learning Model Validation

**Location:** `machine-learning/modelling/polychem_ai.py`

#### Test 8.1: Standalone Script

1. Open terminal in `machine-learning/modelling/`
2. Run:
   ```bash
   python polychem_ai.py
   ```
3. Input SMILES when prompted
4. **Expected:** Prints recommendation without errors

#### Test 8.2: Dataset Loading

1. Check backend startup logs
2. **Expected:** Log shows successful dataset load (row count)

---

## 📋 Critical Issues Log

| Issue                         | Page/Component     | Steps to Reproduce                 | Priority | Status                        |
| ----------------------------- | ------------------ | ---------------------------------- | -------- | ----------------------------- |
| Firebase auth/invalid-api-key | All                | Invalid .env.local keys            | ✅ FIXED | Config fallback added         |
| White screen render           | /                  | Missing Firebase env               | ✅ FIXED | Fallback config + routing fix |
| toFixed crash                 | ChemicalDetailPage | Null molecular_weight              | ✅ FIXED | normalizePrediction guards    |
| Backend 500 on image fail     | /predict API       | Invalid SMILES causing RDKit error | ✅ FIXED | Try-catch with fallback       |
| Missing dependencies          | Backend import     | Missing requests/Pillow            | ✅ FIXED | Added to requirements.txt     |
|                               |                    |                                    |          |                               |

---

## ✅ Final Checklist

- [ ] All **HIGH priority** tests pass (Section 1-2)
- [ ] No console errors in browser DevTools
- [ ] Firebase Console shows user data
- [ ] Backend logs show no uncaught exceptions
- [ ] Images render correctly
- [ ] Dark mode persists after refresh
- [ ] Mobile view looks clean

---

## 🚀 Deployment Readiness

### Before Pushing to Production:

1. **Backend (.env)**
   - [ ] Remove hardcoded API keys from code
   - [ ] Set `GOOGLE_API_KEY` in Koyeb env vars
   - [ ] Verify `/tmp` cache path works on Koyeb

2. **Frontend (.env.local → Vercel)**
   - [ ] Add all `VITE_FIREBASE_*` vars to Vercel env settings
   - [ ] Set `VITE_API_BASE_URL=https://your-backend.koyeb.app`
   - [ ] Test production build: `npm run build && npm run preview`

3. **Database**
   - [ ] Firestore security rules review (currently public?)
   - [ ] Create indexes for `predictions` collection queries

4. **Monitoring**
   - [ ] Setup error logging (Sentry/LogRocket)
   - [ ] Analytics tracking (GA/Posthog)

---

**Notes:**

- Use this checklist before every deployment
- Update **Critical Issues Log** if new bugs found
- Share results with team after each test cycle

**Tested by:** ******\_\_\_******  
**Date:** ******\_\_\_******  
**Overall Result:** ⭕ PASS / ❌ FAIL
