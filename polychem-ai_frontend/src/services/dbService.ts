import { doc, setDoc, updateDoc, arrayUnion, getDoc } from "firebase/firestore"; 
import { db } from "../lib/firebase";
import { type ChemicalData } from "./mockData";

// --- SAVE SEARCH HISTORY ---
export const saveSearchToHistory = async (uid: string, keyword: string) => {
  if (!uid || !keyword.trim()) return;
  try {
    const userRef = doc(db, "users", uid);
    await setDoc(userRef, {
      searchHistory: arrayUnion({
        keyword: keyword,
        timestamp: new Date().toISOString()
      })
    }, { merge: true });
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
    }
    return [];
  } catch (error) {
    console.error("Error getting history:", error); // FIX: Gunakan variabel error
    return [];
  }
};

// --- SAVE TO LIBRARY ---
export const saveToLibrary = async (uid: string, chemical: ChemicalData) => {
  if (!uid || !chemical) return false;
  try {
    const userRef = doc(db, "users", uid);
    await setDoc(userRef, {
      library: arrayUnion({
        ...chemical, 
        savedAt: new Date().toISOString()
      })
    }, { merge: true });
    return true; 
  } catch (error) {
    console.error("Error saving to library:", error);
    return false;
  }
};

// --- REMOVE FROM LIBRARY ---
export const removeFromLibrary = async (uid: string, chemicalId: number) => {
  if (!uid) return false;
  try {
    const userRef = doc(db, "users", uid);
    const docSnap = await getDoc(userRef);

    if (docSnap.exists()) {
      const data = docSnap.data();
      const currentLibrary = data.library || [];

      // Filter item yang akan dihapus
      const newLibrary = currentLibrary.filter((item: ChemicalData) => item.id !== chemicalId);

      await updateDoc(userRef, {
        library: newLibrary
      });
      return true;
    }
    return false;
  } catch (error) {
    console.error("Gagal menghapus dari library:", error); // FIX: Gunakan variabel error
    return false;
  }
};

// --- CHECK IF SAVED ---
export const checkIsSaved = async (uid: string, chemicalId: number) => {
  if (!uid) return false;
  try {
    const userRef = doc(db, "users", uid);
    const docSnap = await getDoc(userRef);
    if (docSnap.exists()) {
      const data = docSnap.data();
      const library = data.library || [];
      return library.some((item: ChemicalData) => item.id === chemicalId);
    }
    return false;
  } catch (error) {
    console.error("Error checking saved status:", error); // FIX: Gunakan variabel error
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
    }
    return [];
  } catch (error) {
    console.error("Gagal mengambil library:", error); // FIX: Gunakan variabel error
    return [];
  }
};