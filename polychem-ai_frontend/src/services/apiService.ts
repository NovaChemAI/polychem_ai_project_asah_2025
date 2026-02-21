// src/services/apiService.ts

const DEFAULT_API_BASE_URL = "https://slim-danika-polychem-ab276767.koyeb.app";
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
    const response = await fetch(`${API_BASE_URL}/predict`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
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
