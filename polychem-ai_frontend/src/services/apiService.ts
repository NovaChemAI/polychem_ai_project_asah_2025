// src/services/apiService.ts

import { auth } from "../lib/firebase";

const DEFAULT_API_BASE_URL =
  "http://localhost:8000";
export const API_BASE_URL = (
  import.meta.env.VITE_API_BASE_URL || DEFAULT_API_BASE_URL
).replace(/\/+$/, "");

interface SimilarCompound {
  name: string;
  smiles: string;
  formula: string;
  molecular_weight: number;
  tg: number;
  image_url: string;
  similarity_percent: number;
}

interface NewCompound {
  name: string;
  smiles: string;
  formula: string;
  molecular_weight: number;
  tg: number;
  tg_justification: string;
  polymer_class: string;
  justifikasi: string;
  image_url: string;
}

export interface PredictionResult {
  status: string;
  input_smiles: string;
  new_compound: NewCompound;
  similar_compounds: SimilarCompound[];
}

const isRecord = (value: unknown): value is Record<string, unknown> => {
  return !!value && typeof value === "object" && !Array.isArray(value);
};

const toStr = (value: unknown, fallback = ""): string => {
  return typeof value === "string" ? value : fallback;
};

const toNum = (value: unknown, fallback = 0): number => {
  const n = typeof value === "number" ? value : Number(value);
  return Number.isFinite(n) ? n : fallback;
};

export const buildAssetUrl = (path: string): string => {
  if (!path) return "";
  if (path.startsWith("http://") || path.startsWith("https://")) return path;
  return `${API_BASE_URL}${path.startsWith("/") ? "" : "/"}${path}`;
};

export const normalizePrediction = (raw: unknown): PredictionResult | null => {
  if (!isRecord(raw)) return null;

  const newRaw = isRecord(raw.new_compound) ? raw.new_compound : {};
  const similarRaw = Array.isArray(raw.similar_compounds)
    ? raw.similar_compounds
    : [];

  const normalized: PredictionResult = {
    status: toStr(raw.status, "success"),
    input_smiles: toStr(raw.input_smiles, toStr(newRaw.smiles)),
    new_compound: {
      name: toStr(newRaw.name, "Unknown Compound"),
      smiles: toStr(newRaw.smiles),
      formula: toStr(newRaw.formula, "-"),
      molecular_weight: toNum(newRaw.molecular_weight, 0),
      tg: toNum(newRaw.tg, 0),
      tg_justification: toStr(
        newRaw.tg_justification,
        "Tidak ada justifikasi Tg.",
      ),
      polymer_class: toStr(newRaw.polymer_class, "Novel Compound"),
      justifikasi: toStr(newRaw.justifikasi, "Tidak ada justifikasi."),
      image_url: toStr(newRaw.image_url),
    },
    similar_compounds: similarRaw.filter(isRecord).map((item) => ({
      name: toStr(item.name, "Unknown Compound"),
      smiles: toStr(item.smiles),
      formula: toStr(item.formula, "-"),
      molecular_weight: toNum(item.molecular_weight, 0),
      tg: toNum(item.tg, 0),
      image_url: toStr(item.image_url),
      similarity_percent: toNum(item.similarity_percent, 0),
    })),
  };

  if (!normalized.new_compound.smiles) {
    return null;
  }

  return normalized;
};

export const predictPolymer = async (
  query: string,
): Promise<PredictionResult | null> => {
  try {
    const headers: Record<string, string> = {
      "Content-Type": "application/json",
    };

    const user = auth.currentUser;
    if (user) {
      const token = await user.getIdToken();
      headers.Authorization = `Bearer ${token}`;
    }

    const response = await fetch(`${API_BASE_URL}/predict`, {
      method: "POST",
      headers,
      body: JSON.stringify({ smiles: query }),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      console.error("Backend Error Detail:", errorData);
      return null;
    }

    const data = await response.json();
    return normalizePrediction(data);
  } catch (error) {
    console.error("Gagal koneksi ke AI:", error);
    return null;
  }
};

// =========================================================
// LIBRARY (via Backend FastAPI, bukan langsung ke Firestore)
// =========================================================

export interface SavedChemicalResponse {
  id: string;
  userId: string;
  name: string;
  smiles: string;
  category: string;
  image?: string;
  image_url?: string;
  properties?: string;
  score?: string | number;
  isAiResult?: boolean;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  savedAt?: any;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
export const saveToLibraryBackend = async (data: any): Promise<boolean> => {
  try {
    const user = auth.currentUser;
    if (!user) {
      console.error("User belum login, tidak ada token.");
      return false;
    }
    const token = await user.getIdToken();

    const response = await fetch(`${API_BASE_URL}/library`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${token}`,
      },
      body: JSON.stringify(data),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      console.error("Gagal simpan ke library:", errorData);
      return false;
    }

    return true;
  } catch (error) {
    console.error("Gagal koneksi saat simpan library:", error);
    return false;
  }
};

export const getUserLibraryBackend = async (): Promise <SavedChemicalResponse[]> => {
  try {
    const user = auth.currentUser;
    if (!user) return [];
    const token = await user.getIdToken();

    const response = await fetch(`${API_BASE_URL}/library`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    });

    if (!response.ok) return [];
    return await response.json();
  } catch (error) {
    console.error("Gagal mengambil library:", error);
    return [];
  }
};

export const removeFromLibraryBackend = async (
  itemId: string,
): Promise<boolean> => {
  try {
    const user = auth.currentUser;
    if (!user) return false;
    const token = await user.getIdToken();

    const response = await fetch(`${API_BASE_URL}/library/${itemId}`, {
      method: "DELETE",
      headers: {
        Authorization: `Bearer ${token}`,
      },
    });

    return response.ok;
  } catch (error) {
    console.error("Gagal menghapus dari library:", error);
    return false;
  }
};

export interface HistoryItemResponse {
  id: string;
  userId: string;
  query: string;
  result_name: string;
  result_smiles: string;
  full_data: string;
  timestamp?: {
    seconds: number;
    nanoseconds?: number;
  };
}

export const getUserHistoryBackend = async (): Promise<HistoryItemResponse[]> => {
  try {
    const user = auth.currentUser;
    if (!user) return [];
    const token = await user.getIdToken();

    const response = await fetch(`${API_BASE_URL}/history/mine`, {
      headers: {
        Authorization: `Bearer ${token}`,
      },
    });

    if (!response.ok) return [];
    return await response.json();
  } catch (error) {
    console.error("Gagal mengambil history:", error);
    return [];
  }
};

export const checkIsSavedBackend = async (
  smiles: string,
  ): Promise<{ isSaved: boolean; itemId: string }> => {
  try {
    const user = auth.currentUser;
    if (!user || !smiles) return { isSaved: false, itemId: "" };
    const token = await user.getIdToken();

    const response = await fetch(
      `${API_BASE_URL}/library/check?smiles=${encodeURIComponent(smiles)}`,
      {
        headers: {
          Authorization: `Bearer ${token}`,
        },
      },
    );

    if (!response.ok) return { isSaved: false, itemId: "" };
    return await response.json();
  } catch (error) {
    console.error("Gagal cek status saved:", error);
    return { isSaved: false, itemId: "" };
  }
};


export const syncUserProfile = async (
  name?: string,
  photoUrl?: string,
): Promise<boolean> => {
  try {
    const user = auth.currentUser;
    if (!user) return false;
    const token = await user.getIdToken();

    const response = await fetch(`${API_BASE_URL}/auth/sync-profile`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${token}`,
      },
      body: JSON.stringify({ name, photo_url: photoUrl }),
    });

    return response.ok;
  } catch (error) {
    console.error("Gagal sinkronisasi profil:", error);
    return false;
  }
};


// =========================================================
// CAPTCHA VERIFICATION
// =========================================================

export const verifyCaptcha = async (captchaToken: string): Promise<boolean> => {
  try {
    const response = await fetch(`${API_BASE_URL}/auth/verify-captcha`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ captcha_token: captchaToken }),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      console.error("Verifikasi captcha gagal:", errorData);
      return false;
    }

    return true;
  } catch (error) {
    console.error("Gagal koneksi saat verifikasi captcha:", error);
    return false;
  }
};