import { db } from '../lib/firebase';
import { 
  doc, 
  setDoc, 
  deleteDoc, 
  getDoc, 
  collection, 
  getDocs, 
  query, 
  where,
  serverTimestamp,
  addDoc,
  orderBy,
  limit
} from 'firebase/firestore';

// Nama Koleksi
const COLLECTION_NAME = "saved_chemicals";
const HISTORY_COLLECTION = "search_history";

// --- INTERFACES ---
export interface SavedChemical {
  id: string;
  userId: string;
  name: string;
  smiles: string;
  category: string;
  
  // PERBAIKAN: Support kedua nama properti
  image?: string;       // Legacy
  image_url?: string;   // Sesuai Backend Python

  properties?: string; 
  score?: string | number;
  isAiResult?: boolean;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  savedAt?: any; 
}

export interface HistoryItem {
  id: string;
  userId: string;
  query: string;
  result_name: string;
  result_smiles: string;
  full_data: string;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  timestamp: any;
}

// --- LIBRARY FUNCTIONS (Manual Save) ---

// eslint-disable-next-line @typescript-eslint/no-explicit-any
export const saveToLibrary = async (userId: string, data: any): Promise<boolean> => {
  try {
    let docId: string;
    if (data.id && !data.isAiResult) {
        docId = data.id.toString(); 
    } else {
        docId = `ai_${Date.now()}_${Math.floor(Math.random() * 1000)}`;
    }

    let cleanProperties = "";
    if (data.isAiResult && typeof data.properties === 'object') {
        cleanProperties = JSON.stringify(data.properties);
    } else if (typeof data.properties === 'string') {
        cleanProperties = data.properties;
    } else {
        cleanProperties = JSON.stringify(data.properties || {});
    }

    // --- PERBAIKAN PENTING DI SINI ---
    // Backend Python mengirim 'image_url'. Kita harus menangkapnya.
    // Kita simpan ke field 'image_url' di Firestore juga agar konsisten.
    const imageUrlToSave = data.image_url || data.image || "";

    const payload: SavedChemical = {
        id: docId,
        userId: userId,
        name: data.name || "Unknown Compound",
        smiles: data.smiles || "",
        category: data.category || "Uncategorized",
        
        // Simpan data gambar
        image: imageUrlToSave,      // Backward compatibility (opsional)
        image_url: imageUrlToSave,  // Field utama untuk masa depan
        
        properties: cleanProperties,
        score: data.score || "N/A",
        isAiResult: true,
        savedAt: serverTimestamp()
    };

    const compositeId = `${userId}_${docId}`;
    await setDoc(doc(db, COLLECTION_NAME, compositeId), payload);
    return true;
  } catch (error) {
    console.error("Error saving to library:", error);
    return false;
  }
};

export const removeFromLibrary = async (userId: string, itemId: string | number): Promise<boolean> => {
  try {
    const docId = itemId.toString();
    const compositeId = `${userId}_${docId}`;
    await deleteDoc(doc(db, COLLECTION_NAME, compositeId));
    return true;
  } catch (error) {
    console.error("Error removing:", error);
    return false;
  }
};

export const checkIsSaved = async (userId: string, itemId: string | number): Promise<boolean> => {
  try {
    if (!itemId) return false;
    const docId = itemId.toString();
    const compositeId = `${userId}_${docId}`;
    const docRef = doc(db, COLLECTION_NAME, compositeId);
    const docSnap = await getDoc(docRef);
    return docSnap.exists();
  } catch (error) {
    console.error("Error checking status:", error);
    return false;
  }
};

export const getUserLibrary = async (userId: string): Promise<SavedChemical[]> => {
  try {
    const q = query(collection(db, COLLECTION_NAME), where("userId", "==", userId));
    const querySnapshot = await getDocs(q);
    const items: SavedChemical[] = [];
    querySnapshot.forEach((doc) => items.push(doc.data() as SavedChemical));
    return items;
  } catch (error) {
    console.error("Error getting library:", error);
    return [];
  }
};

// --- HISTORY FUNCTIONS (Automatic) ---

// eslint-disable-next-line @typescript-eslint/no-explicit-any
export const addToHistory = async (userId: string, queryText: string, resultData: any) => {
  try {
    console.log("Mencoba menyimpan history..."); 
    const nc = resultData.new_compound || {};
    
    const payload = {
        userId: userId,
        query: queryText,
        result_name: nc.name || "Unknown",
        result_smiles: nc.smiles || "",
        full_data: JSON.stringify(resultData),
        timestamp: serverTimestamp()
    };

    await addDoc(collection(db, HISTORY_COLLECTION), payload);
    console.log("✅ History berhasil disimpan ke Firebase");
  } catch (error) {
    console.error("❌ Error logging history:", error);
  }
};

export const getUserHistory = async (userId: string): Promise<HistoryItem[]> => {
  try {
    const q = query(
        collection(db, HISTORY_COLLECTION), 
        where("userId", "==", userId),
        orderBy("timestamp", "desc"),
        limit(50)
    );
    
    const querySnapshot = await getDocs(q);
    const items: HistoryItem[] = [];
    
    querySnapshot.forEach((doc) => {
      const d = doc.data();
      items.push({
        id: doc.id,
        userId: d.userId,
        query: d.query,
        result_name: d.result_name,
        result_smiles: d.result_smiles,
        full_data: d.full_data,
        timestamp: d.timestamp
      });
    });

    return items;
  } catch (error) {
    console.error("Error getting history:", error);
    return [];
  }
};