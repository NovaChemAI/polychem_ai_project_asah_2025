# 🚀 Production Deployment Guide - Polychem AI

**Last Updated:** February 22, 2026  
**Deployment Stack:** Backend → Koyeb | Frontend → Vercel | Auth/DB → Firebase

---

## 📋 Pre-Deployment Checklist

### ✅ Code Quality & Security

- [x] No hardcoded API keys in code
- [x] `.env.local` in `.gitignore`
- [x] CORS configured with environment variable
- [x] Input validation (Pydantic max_length=500)
- [x] Error sanitization (no internal details leaked)
- [x] Firebase fallback config for dev
- [x] XSS protection (React auto-escaping, no dangerouslySetInnerHTML)
- [x] Cache directory writable in container (/tmp)
- [x] Static files path Koyeb-safe (/tmp/static)

### ✅ Performance

- [x] LLM caching (diskcache)
- [x] Image generation caching (MD5 hash-based)
- [x] Dataset loaded once at startup
- [x] Response normalization (frontend guards null)
- [x] Lazy loading untuk images

### ⚠️ Security Improvements Needed (Optional)

- [ ] Rate limiting (see Recommendations section)
- [ ] Firestore security rules (see Firebase section)
- [ ] API key header authentication (if public access concern)
- [ ] Request size limits (already at 500 char SMILES)
- [ ] HTTPS enforcement (handled by Koyeb/Vercel)

---

## 🔧 Backend Deployment (Koyeb)

### 1. Prerequisites

- **Koyeb Account:** https://koyeb.com (free tier available)
- **Google Gemini API Key:** https://aistudio.google.com/apikey
- **Git Repository:** GitHub/GitLab repo with backend code

### 2. Koyeb Setup

#### A. Create New Service

1. Login to Koyeb → **Create Service**
2. **Deployment Method:** GitHub (connect your repo)
3. **Branch:** `main` or `master`
4. **Builder:** Dockerfile
5. **Dockerfile Path:** `polychem-ai_backend/Dockerfile`
6. **Build Context:** `polychem-ai_backend/`

#### B. Configure Service Settings

**Instance Type:**

- **Size:** `Nano` (512MB RAM) for testing, `Small` (1GB) for production
- **Regions:** Choose closest to your users (e.g., `fra` for Europe, `was` for US)
- **Scaling:** 1 instance (free tier) or enable auto-scaling

**Ports:**

- **Port:** `8000` (auto-detected from Dockerfile)
- **Protocol:** HTTP
- **Public:** ✅ Enabled

#### C. Environment Variables (Critical!)

**In Koyeb Dashboard → Service → Settings → Environment Variables:**

| Variable         | Value                         | Required       |
| ---------------- | ----------------------------- | -------------- |
| `GOOGLE_API_KEY` | Your Gemini API key           | ✅ **YES**     |
| `FRONTEND_URL`   | `https://your-app.vercel.app` | ⚠️ Recommended |
| `CACHE_VERSION`  | `v1`                          | Optional       |
| `HISTORY_LIMIT`  | `50`                          | Optional       |

**IMPORTANT:** Never commit `GOOGLE_API_KEY` to git!

#### D. Deploy

1. Click **Deploy**
2. Wait 3-5 minutes for build
3. Check logs for `Dataset loaded!`
4. Copy your backend URL: `https://your-app.koyeb.app`

### 3. Verify Backend

```bash
# Health check
curl https://your-app.koyeb.app/health
# Expected: {"status":"ok"}

# Test prediction
curl -X POST https://your-app.koyeb.app/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO"}'
# Expected: JSON with new_compound + similar_compounds
```

### 4. Koyeb-Specific Notes

**✅ What's Already Handled:**

- `/tmp` directory writable (cache + images persisted during container lifetime)
- `PORT` env injected automatically (Dockerfile uses `${PORT:-8000}`)
- HTTPS auto-configured with Let's Encrypt
- `--proxy-headers` enabled for correct client IP forwarding

**⚠️ Limitations:**

- Cache resets on container restart (acceptable for LLM responses)
- `/app` directory is read-only (don't write files there)
- Free tier: 1 instance max, container sleep after 15min inactivity

---

## 🎨 Frontend Deployment (Vercel)

### 1. Prerequisites

- **Vercel Account:** https://vercel.com (free tier generous)
- **Firebase Project:** https://console.firebase.google.com
- **Backend URL:** From Koyeb deployment above

### 2. Firebase Setup

#### A. Create Firebase Project

1. Go to Firebase Console → **Add Project**
2. **Name:** `polychem-ai` (or your choice)
3. **Analytics:** Optional (can disable)

#### B. Enable Authentication

1. **Build → Authentication → Get Started**
2. **Sign-in method:**
   - ✅ **Email/Password** (enable)
   - ✅ **Google** (enable, configure OAuth consent)
3. **Authorized domains:** Add your Vercel domain when ready

#### C. Enable Firestore

1. **Build → Firestore Database → Create database**
2. **Mode:** Start in **production mode** (we'll add rules next)
3. **Location:** Choose closest to users

#### D. Security Rules (IMPORTANT!)

**Replace default Firestore rules with:**

```javascript
rules_version = '2';
service cloud.firestore {
  match /databases/{database}/documents {

    // Users can only read/write their own predictions
    match /predictions/{predictionId} {
      allow read: if request.auth != null && request.auth.uid == resource.data.userId;
      allow create: if request.auth != null && request.auth.uid == request.resource.data.userId;
      allow update, delete: if request.auth != null && request.auth.uid == resource.data.userId;
    }

    // Users collection (optional, for profiles)
    match /users/{userId} {
      allow read, write: if request.auth != null && request.auth.uid == userId;
    }

    // Deny all other access
    match /{document=**} {
      allow read, write: if false;
    }
  }
}
```

**Deploy rules:** Click **Publish** in Firestore Rules editor

#### E. Get Firebase Config

1. **Project Settings** (gear icon) → **General**
2. **Your apps** → **Web app** (</> icon)
3. **Register app:** Name it `Polychem AI Web`
4. Copy **firebaseConfig** object values

### 3. Vercel Deployment

#### A. Connect Repository

1. Login to Vercel → **Add New Project**
2. **Import Git Repository:** Select your repo
3. **Root Directory:** `polychem-ai_frontend`
4. **Framework Preset:** Vite

#### B. Environment Variables

**In Vercel → Project Settings → Environment Variables:**

| Variable                            | Value                          | Source            |
| ----------------------------------- | ------------------------------ | ----------------- |
| `VITE_FIREBASE_API_KEY`             | From Firebase config           | Firebase Console  |
| `VITE_FIREBASE_AUTH_DOMAIN`         | `your-project.firebaseapp.com` | Firebase Console  |
| `VITE_FIREBASE_PROJECT_ID`          | `your-project`                 | Firebase Console  |
| `VITE_FIREBASE_STORAGE_BUCKET`      | `your-project.appspot.com`     | Firebase Console  |
| `VITE_FIREBASE_MESSAGING_SENDER_ID` | `123456789`                    | Firebase Console  |
| `VITE_FIREBASE_APP_ID`              | `1:123:web:abc123`             | Firebase Console  |
| `VITE_API_BASE_URL`                 | `https://your-app.koyeb.app`   | Koyeb backend URL |

**Environment Scope:** All (Production, Preview, Development)

#### C. Build Settings

**Should auto-detect, but verify:**

- **Build Command:** `npm run build`
- **Output Directory:** `dist`
- **Install Command:** `npm install`
- **Development Command:** `npm run dev`

#### D. Deploy

1. Click **Deploy**
2. Wait 2-3 minutes
3. Vercel assigns URL: `https://your-app.vercel.app`

### 4. Update Backend CORS

**In Koyeb Environment Variables:**

- Add/Update `FRONTEND_URL` = `https://your-app.vercel.app`
- Redeploy backend to apply changes

### 5. Verify Frontend

1. Open `https://your-app.vercel.app`
2. **Test Registration:**
   - Click **Register**
   - Fill form → Submit
   - Check Firebase Console → Authentication (user should appear)
3. **Test Prediction:**
   - Input SMILES: `CCO`
   - Check network tab: Should hit `https://your-app.koyeb.app/predict`
   - Result page should load with images

---

## 🔐 Security Hardening (Post-Deployment)

### 1. Backend Rate Limiting (Recommended)

**Add to `requirements.txt`:**

```
slowapi==0.1.9
```

**Add to `app/main.py`:**

```python
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

limiter = Limiter(key_func=get_remote_address, default_limits=["60/minute"])
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

@app.post("/predict")
@limiter.limit("10/minute")  # Max 10 predictions per minute per IP
def predict(req: RecommendRequest, request: Request):
    # ... existing code
```

**Why:** Prevent abuse, protect LLM API quota

### 2. Firestore Indexes (Performance)

**In Firestore Console → Indexes:**

Create composite index:

- **Collection:** `predictions`
- **Fields:** `userId` (Ascending), `timestamp` (Descending)

**Why:** Speeds up "My History" queries

### 3. Firebase Authentication Restrictions

**In Firebase Console → Authentication → Settings:**

**Authorized Domains:**

- Remove `localhost` (only for dev)
- Keep only your Vercel domain

**Email Enumeration Protection:**

- ✅ Enable (prevents email harvesting)

### 4. Content Security Policy (Vercel)

**Create `polychem-ai_frontend/vercel.json`:**

```json
{
  "rewrites": [{ "source": "/(.*)", "destination": "/index.html" }],
  "headers": [
    {
      "source": "/(.*)",
      "headers": [
        {
          "key": "X-Content-Type-Options",
          "value": "nosniff"
        },
        {
          "key": "X-Frame-Options",
          "value": "DENY"
        },
        {
          "key": "X-XSS-Protection",
          "value": "1; mode=block"
        },
        {
          "key": "Referrer-Policy",
          "value": "strict-origin-when-cross-origin"
        }
      ]
    }
  ]
}
```

**Redeploy frontend to apply**

---

## 📊 Monitoring & Maintenance

### Backend (Koyeb)

**Logs:**

- Dashboard → Service → Logs
- Check for `Dataset loaded!` on startup
- Monitor LLM errors (rate limits, timeouts)

**Metrics:**

- CPU/Memory usage (upgrade if >80% sustained)
- Request count
- Response times

**Rebuild Triggers:**

- Update Python dependencies
- Change environment variables
- New dataset

### Frontend (Vercel)

**Analytics:**

- Vercel → Analytics (Web Vitals)
- Monitor Core Web Vitals (LCP, FID, CLS)

**Error Tracking (Optional):**

- Integrate Sentry: `npm install @sentry/react`
- Add to `main.tsx`

**Logs:**

- Vercel → Deployments → Functions (for SSR if any)
- Browser Console for client errors

### Firebase

**Usage:**

- Console → Usage and billing
- Monitor authentication quota
- Firestore reads/writes

**Security:**

- Review authentication logs
- Check for suspicious activity

---

## 🚨 Common Deployment Issues

### Backend Issues

**❌ "GOOGLE_API_KEY belum diset"**

- **Fix:** Add `GOOGLE_API_KEY` in Koyeb env vars
- Redeploy service

**❌ CORS error from frontend**

- **Fix:** Set `FRONTEND_URL` in Koyeb env to match Vercel URL
- Ensure frontend uses correct API base URL

**❌ "Permission denied" writing cache**

- **Fix:** Check `CACHE_DIR` and `STATIC_DIR` use `/tmp/` prefix
- Already configured correctly in code

**❌ Port binding error**

- **Fix:** Dockerfile uses `${PORT:-8000}` (auto-handled by Koyeb)

### Frontend Issues

**❌ Firebase "invalid-api-key"**

- **Fix:** Verify all `VITE_FIREBASE_*` env vars in Vercel
- Redeploy after adding vars

**❌ "Failed to fetch" from backend**

- **Fix:** Check `VITE_API_BASE_URL` points to correct Koyeb URL
- Verify backend `/health` endpoint returns 200

**❌ Blank page (white screen)**

- **Fix:** Check browser console for errors
- Verify routing config in `App.tsx`
- Check Firebase config fallback

**❌ Images not loading**

- **Fix:** Backend must be running
- Check CORS headers on `/static/` path
- Verify `buildAssetUrl()` constructs correct URL

---

## 🎯 Production Optimization

### Backend

**1. Dataset Caching:**

```python
# Already implemented in store.py
# Dataset loads once at startup (lifespan event)
```

**2. LLM Response Caching:**

```python
# Already implemented in llm.py + cache.py
# Uses diskcache with TTL
```

**3. Image Caching:**

```python
# Already implemented in images.py
# MD5 hash prevents regeneration
```

### Frontend

**1. Code Splitting:**

```typescript
// Use React.lazy for route-based splitting
const LibraryPage = lazy(() => import("./pages/LibraryPage"));
```

**2. Image Optimization:**

```typescript
// Add loading="lazy" to <img> tags in ChemicalDetailPage
<img src={buildAssetUrl(compound.image_url)} loading="lazy" />
```

**3. Bundle Size:**

```bash
# Analyze bundle in Vercel deployment logs
# Consider removing unused dependencies
```

---

## 📝 Deployment Commands Reference

### Backend (Local Testing)

```bash
# Install dependencies
cd polychem-ai_backend
pip install -r requirements.txt

# Run locally
export GOOGLE_API_KEY=your_key_here
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

### Frontend (Local Testing)

```bash
# Install dependencies
cd polychem-ai_frontend
npm install

# Create .env.local (copy from .env.example)
cp .env.example .env.local
# Edit .env.local with Firebase keys

# Run dev server
npm run dev

# Production build test
npm run build
npm run preview
```

### Git Workflow

```bash
# Before pushing to production
git status
# Ensure .env, .env.local not staged

git add .
git commit -m "Production ready: security + deployment config"
git push origin main
# Koyeb + Vercel auto-deploy on push
```

---

## ✅ Final Verification Checklist

**Before Going Live:**

### Security

- [ ] No API keys in git history
- [ ] CORS restricted to production domain
- [ ] Firestore security rules deployed
- [ ] Firebase authorized domains configured
- [ ] `.env.local` in `.gitignore`

### Functionality

- [ ] Backend `/health` returns 200
- [ ] Backend `/predict` accepts valid SMILES
- [ ] Frontend loads without errors
- [ ] User registration/login works
- [ ] SMILES prediction returns results
- [ ] Images load correctly
- [ ] Similar compounds display
- [ ] Library save/retrieve works
- [ ] History shows past predictions

### Performance

- [ ] LCP < 2.5s (Vercel Analytics)
- [ ] Backend response < 5s for predictions
- [ ] Images cached (check DevTools Network)
- [ ] Dataset loads once at startup

### Monitoring

- [ ] Koyeb logs visible
- [ ] Vercel deployment successful
- [ ] Firebase usage tracking enabled
- [ ] Error tracking configured (optional)

---

## 🆘 Support & Resources

**Koyeb Docs:** https://koyeb.com/docs  
**Vercel Docs:** https://vercel.com/docs  
**Firebase Docs:** https://firebase.google.com/docs  
**FastAPI Docs:** https://fastapi.tiangolo.com  
**React Docs:** https://react.dev

**Common Issues:**

- Check [MANUAL_TESTING_CHECKLIST.md](MANUAL_TESTING_CHECKLIST.md) for testing procedures
- Review [BACKEND_COMPREHENSIVE_REPORT.md](BACKEND_COMPREHENSIVE_REPORT.md) for architecture
- See [TESTING_GUIDE.md](TESTING_GUIDE.md) for validation

---

**🚀 Deployment Complete! Your Polychem AI is production-ready.**

**Next Steps:**

1. Monitor usage for 48 hours
2. Implement rate limiting if abuse detected
3. Set up error tracking (Sentry)
4. Add Google Analytics (optional)
5. Create user documentation
