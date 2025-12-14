// src/services/apiService.ts

// URL Backend (Port 8000 sesuai terminal Anda)
const API_BASE_URL = "http://127.0.0.1:8000";

export interface PredictionResult {
  // Kita definisikan field yang kemungkinan besar ada (opsional)
  smiles?: string;
  biodegradability_score?: number;

  // Baris sakti untuk mematikan error 'any' khusus di sini:
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  [key: string]: any; 
}

export const predictPolymer = async (query: string): Promise<PredictionResult | null> => {
  try {
    // Kita coba nembak ke endpoint /predict
    const response = await fetch(`${API_BASE_URL}/predict`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      // Backend FastAPI biasanya minta body JSON. Kita kirim key 'smiles'.
      body: JSON.stringify({ smiles: query }), 
    });

    if (!response.ok) {
      // Jika error, kita baca pesan errornya dari backend
      const errorData = await response.json().catch(() => ({}));
      console.error("Backend Error Detail:", errorData);
      throw new Error(`API Error: ${response.status}`);
    }

    const data = await response.json();
    return data;

  } catch (error) {
    console.error("Gagal koneksi ke AI:", error);
    return null;
  }
};
