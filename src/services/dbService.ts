// Ganti import: Hapus 'updateDoc', Tambahkan 'setDoc'
import { doc, setDoc, arrayUnion, getDoc } from "firebase/firestore"; 
import { db } from "../lib/firebase";
import { type ChemicalData } from "./mockData";

// --- SAVE SEARCH HISTORY ---
export const saveSearchToHistory = async (uid: string, keyword: string) => {
  if (!uid || !keyword.trim()) return;

  try {
    const userRef = doc(db, "users", uid);
    
    // GUNAKAN setDoc + merge: true (Lebih Aman daripada updateDoc)
    await setDoc(userRef, {
      searchHistory: arrayUnion({
        keyword: keyword,
        timestamp: new Date().toISOString()
      })
    }, { merge: true }); // <--- Opsi merge: true sangat penting!
    
    console.log(`History saved: ${keyword}`);
  } catch (error) {
    console.error("Error saving history:", error);
  }
};

// --- GET USER HISTORY ---
export const getUserHistory = async (uid: string) => {
  if (!uid) return [];

  try {
    const userRef = doc(db, "users", uid);
    const docSnap = await getDoc(userRef);

    if (docSnap.exists()) {
      const data = docSnap.data();
      return data.searchHistory ? data.searchHistory.reverse() : [];
    } else {
      return [];
    }
  } catch (error) {
    console.error("Gagal mengambil history:", error);
    return [];
  }
};

// --- SAVE TO LIBRARY (LOGIC FIX) ---
export const saveToLibrary = async (uid: string, chemical: ChemicalData) => {
  if (!uid || !chemical) return false;

  try {
    const userRef = doc(db, "users", uid);
    
    // FIX: Ganti updateDoc menjadi setDoc dengan merge: true
    // Ini akan otomatis MEMBUAT dokumen user jika belum ada.
    await setDoc(userRef, {
      library: arrayUnion({
        ...chemical, 
        savedAt: new Date().toISOString()
      })
    }, { merge: true });
    
    console.log(`Chemical saved to Library: ${chemical.name}`);
    return true; 
  } catch (error) {
    console.error("Error saving to library:", error);
    return false;
  }
};

// --- GET LIBRARY ---
export const getLibrary = async (uid: string) => {
  if (!uid) return [];

  try {
    const userRef = doc(db, "users", uid);
    const docSnap = await getDoc(userRef);

    if (docSnap.exists()) {
      const data = docSnap.data();
      return data.library ? data.library.reverse() : [];
    } else {
      return [];
    }
  } catch (error) {
    console.error("Gagal mengambil library:", error);
    return [];
  }
};