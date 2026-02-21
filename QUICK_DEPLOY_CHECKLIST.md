# 🚀 Quick Deploy Checklist - Polychem AI

**Use this for fast deployment. For detailed guide, see [PRODUCTION_DEPLOYMENT_GUIDE.md](PRODUCTION_DEPLOYMENT_GUIDE.md)**

---

## 📋 Pre-Deployment (5 minutes)

### 1. Get API Keys

- [ ] **Google Gemini API Key**
  - Go to: https://aistudio.google.com/apikey
  - Click "Create API Key"
  - Copy key (starts with `AIzaSy...`)

- [ ] **Firebase Project**
  - Go to: https://console.firebase.google.com
  - Create new project: `polychem-ai`
  - Enable Authentication (Email + Google)
  - Enable Firestore Database
  - Get config from Project Settings → Web App

---

## 🔧 Backend Deployment (Koyeb) - 10 minutes

### 1. Create Koyeb Account

- Sign up: https://koyeb.com (free tier)
- Connect GitHub account

### 2. Deploy Backend

1. **Create Service** → **Deploy from GitHub**
2. **Repository:** Select your repo
3. **Builder:** Dockerfile
4. **Dockerfile path:** `polychem-ai_backend/Dockerfile`
5. **Build context:** `polychem-ai_backend/`

### 3. Environment Variables

Click **Add Environment Variable**:

```bash
GOOGLE_API_KEY=<paste_your_gemini_key_here>
FRONTEND_URL=https://your-app.vercel.app  # Will add after frontend deploy
```

### 4. Deploy & Get URL

- Click **Deploy** (wait 3-5 minutes)
- Copy backend URL: `https://your-app-name.koyeb.app`
- Test: `curl https://your-app-name.koyeb.app/health`
  - Should return: `{"status":"ok"}`

---

## 🎨 Frontend Deployment (Vercel) - 10 minutes

### 1. Create Vercel Account

- Sign up: https://vercel.com
- Connect GitHub account

### 2. Deploy Frontend

1. **Add New Project** → **Import Git Repository**
2. **Root Directory:** `polychem-ai_frontend`
3. **Framework Preset:** Vite (auto-detected)

### 3. Environment Variables

Add these (**CRITICAL - copy from Firebase Console**):

```bash
# From Firebase Console → Project Settings → General → Your apps → Web app
VITE_FIREBASE_API_KEY=AIzaSy...
VITE_FIREBASE_AUTH_DOMAIN=your-project.firebaseapp.com
VITE_FIREBASE_PROJECT_ID=your-project
VITE_FIREBASE_STORAGE_BUCKET=your-project.appspot.com
VITE_FIREBASE_MESSAGING_SENDER_ID=123456789
VITE_FIREBASE_APP_ID=1:123:web:abc123

# Your Koyeb backend URL from step above
VITE_API_BASE_URL=https://your-app-name.koyeb.app
```

### 4. Deploy & Get URL

- Click **Deploy** (wait 2-3 minutes)
- Copy frontend URL: `https://your-app.vercel.app`

---

## 🔐 Firebase Configuration - 5 minutes

### 1. Deploy Firestore Security Rules

1. Firebase Console → **Firestore Database** → **Rules** tab
2. Copy content from `firestore.rules` file in project root
3. Click **Publish**

### 2. Configure Authorized Domains

1. Firebase Console → **Authentication** → **Settings** → **Authorized domains**
2. Click **Add domain**
3. Add: `your-app.vercel.app` (your Vercel URL)

### 3. Update Backend CORS

1. Go to Koyeb → Your Service → **Environment**
2. Update `FRONTEND_URL` variable: `https://your-app.vercel.app`
3. Click **Redeploy**

---

## ✅ Verification - 5 minutes

### 1. Test Backend

```bash
# Health check
curl https://your-app-name.koyeb.app/health

# Prediction test
curl -X POST https://your-app-name.koyeb.app/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO"}'
```

Expected: JSON response with `new_compound` and `similar_compounds`

### 2. Test Frontend

1. Open: `https://your-app.vercel.app`
2. **Register** new account (test email)
3. **Login** with test account
4. **Predict** with SMILES: `CCO` or `c1ccccc1`
5. Verify:
   - [ ] Results load
   - [ ] Images display
   - [ ] Similar compounds show
   - [ ] Can save to library

### 3. Check Browser Console

- Open DevTools (F12)
- Should have NO errors (Firebase warnings OK)

---

## 🚨 Troubleshooting

### Backend Issues

**❌ "GOOGLE_API_KEY belum diset"**

- Fix: Add `GOOGLE_API_KEY` in Koyeb environment variables
- Redeploy service

**❌ Build fails**

- Check Koyeb logs for errors
- Verify Dockerfile path is correct
- Ensure all files committed to git

### Frontend Issues

**❌ Firebase "invalid-api-key"**

- Verify all `VITE_FIREBASE_*` vars in Vercel
- Check for typos in API keys
- Redeploy after adding vars

**❌ CORS error**

- Ensure `FRONTEND_URL` in Koyeb matches Vercel URL exactly
- Include `https://` in URL
- Redeploy backend after changing

**❌ "Failed to fetch"**

- Check `VITE_API_BASE_URL` in Vercel env
- Verify backend is running (test /health endpoint)
- Check backend logs for errors

---

## 📊 Post-Deployment

### Monitor (First 24 Hours)

**Koyeb Dashboard:**

- Check logs for errors
- Monitor CPU/Memory usage
- Verify requests coming through

**Vercel Dashboard:**

- Check Analytics → Web Vitals
- Monitor deployment logs
- Check function invocations

**Firebase Console:**

- Verify users can register/login
- Check Firestore for saved predictions
- Monitor usage and billing

---

## 🎯 Quick Reference

| Service      | URL                               | Dashboard                           |
| ------------ | --------------------------------- | ----------------------------------- |
| **Backend**  | `https://your-app-name.koyeb.app` | https://app.koyeb.com               |
| **Frontend** | `https://your-app.vercel.app`     | https://vercel.com/dashboard        |
| **Firebase** | N/A                               | https://console.firebase.google.com |

### Important Files

- **Backend Config:** `polychem-ai_backend/.env.example`
- **Frontend Config:** `polychem-ai_frontend/.env.example`
- **Security Rules:** `firestore.rules`
- **Deployment Guide:** `PRODUCTION_DEPLOYMENT_GUIDE.md`

---

## ✅ All Done!

**Total Time:** ~35 minutes  
**Status:** Production Ready 🚀

**Next Steps:**

1. Share app link with users
2. Monitor usage for first week
3. Add rate limiting if needed (see full guide)
4. Set up error tracking (optional: Sentry)

**Need Help?** See [PRODUCTION_DEPLOYMENT_GUIDE.md](PRODUCTION_DEPLOYMENT_GUIDE.md) for detailed troubleshooting.
