import { initializeApp } from "firebase/app";
import { getAuth, GoogleAuthProvider } from "firebase/auth";
import { getFirestore } from "firebase/firestore";

const fallbackApiKey = "AIzaSyDUMMYKEYDUMMYKEYDUMMYKEYDUMMYKEY";

const firebaseConfig = {
  apiKey: import.meta.env.VITE_FIREBASE_API_KEY || fallbackApiKey,
  authDomain:
    import.meta.env.VITE_FIREBASE_AUTH_DOMAIN || "local-dev.firebaseapp.com",
  projectId: import.meta.env.VITE_FIREBASE_PROJECT_ID || "local-dev",
  storageBucket:
    import.meta.env.VITE_FIREBASE_STORAGE_BUCKET || "local-dev.appspot.com",
  messagingSenderId:
    import.meta.env.VITE_FIREBASE_MESSAGING_SENDER_ID || "000000000000",
  appId: import.meta.env.VITE_FIREBASE_APP_ID || "1:000000000000:web:localdev",
};

if (!import.meta.env.VITE_FIREBASE_API_KEY) {
  console.warn(
    "Firebase env belum diset. App tetap jalan dengan local fallback config.",
  );
}

const app = initializeApp(firebaseConfig);

export const auth = getAuth(app);
export const googleProvider = new GoogleAuthProvider();
export const db = getFirestore(app);
