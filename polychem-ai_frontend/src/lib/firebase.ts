// src/lib/firebase.ts
import { initializeApp } from "firebase/app";
import { getAuth, GoogleAuthProvider } from "firebase/auth";
import { getFirestore } from "firebase/firestore";

const firebaseConfig = {
  apiKey: "AIzaSyDjmkdAmsHWi8CD7R_OE6yhK0E1kpH68RI",
  authDomain: "novachemai.firebaseapp.com",
  projectId: "novachemai",
  storageBucket: "novachemai.firebasestorage.app",
  messagingSenderId: "423081510336",
  appId: "1:423081510336:web:7a29e004711b7eb836b425"
};
// Selesai Bagian Config 

// 1. Inisialisasi Aplikasi
const app = initializeApp(firebaseConfig);

// 2. Siapkan Fitur Login (Auth)
export const auth = getAuth(app);
export const googleProvider = new GoogleAuthProvider();

// 3. Siapkan Database (Firestore)
export const db = getFirestore(app);