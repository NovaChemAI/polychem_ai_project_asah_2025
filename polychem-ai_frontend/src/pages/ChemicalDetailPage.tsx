import { useLocation, useNavigate } from 'react-router-dom';
import { useEffect, useState } from 'react';
import toast from 'react-hot-toast';
import { saveToLibrary } from '../services/dbService';
import { auth } from '../lib/firebase';

// 1. DEFINISI TIPE DATA (Interface)
interface SimilarCompound {
  name: string;
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

interface PredictionData {
  new_compound: NewCompound;
  similar_compounds: SimilarCompound[];
}

function ChemicalDetailPage() {
  const location = useLocation();
  const navigate = useNavigate();
  const [isSaving, setIsSaving] = useState(false);

  // 2. GUNAKAN TIPE DATA DI SINI
  const predictionData = location.state?.predictionData as PredictionData;
  
  // Safe access (cegah crash jika data null)
  const compound = predictionData?.new_compound;
  const similar = predictionData?.similar_compounds || [];

  useEffect(() => {
    if (!predictionData) {
      navigate('/');
    }
  }, [predictionData, navigate]);

  if (!predictionData) return null;

  const handleSave = async () => {
    const user = auth.currentUser;
    if (!user) {
      toast.error("Silakan login untuk menyimpan.");
      return;
    }
    setIsSaving(true);
    
    const payload = {
      ...compound,
      isAiResult: true,
      category: compound.polymer_class || "Novel Compound"
    };

    const success = await saveToLibrary(user.uid, payload);
    if (success) toast.success("Berhasil disimpan ke Library!");
    else toast.error("Gagal menyimpan.");
    
    setIsSaving(false);
  };

  return (
    <div className="max-w-6xl mx-auto px-4 py-8">
      <button onClick={() => navigate(-1)} className="mb-6 text-gray-500 hover:text-blue-600 flex items-center gap-2 transition-colors">
        ← Kembali
      </button>

      {/* --- KARTU UTAMA --- */}
      <div className="bg-card rounded-2xl shadow-xl border border-border overflow-hidden flex flex-col md:flex-row">
        <div className="md:w-1/3 bg-gray-50 dark:bg-slate-800 p-8 flex items-center justify-center border-r border-border">
          {compound.image_url ? (
            <img 
              src={`http://127.0.0.1:8000${compound.image_url}`} 
              alt="Structure" 
              className="w-full h-auto max-w-[300px] object-contain mix-blend-multiply dark:mix-blend-normal filter dark:invert"
              onError={(e) => {
                (e.target as HTMLImageElement).src = 'https://via.placeholder.com/300?text=No+Structure';
              }}
            />
          ) : (
            <div className="text-gray-400">No Image</div>
          )}
        </div>

        <div className="p-8 md:w-2/3 flex flex-col">
          <div className="flex justify-between items-start mb-4">
            <div>
              <span className="inline-block px-3 py-1 text-xs font-bold tracking-wider text-blue-600 uppercase bg-blue-100 rounded-full mb-2">
                {compound.polymer_class || "Novel Compound"}
              </span>
              <h1 className="text-3xl md:text-4xl font-extrabold text-main mb-2">
                {/* Logic Nama: Backend Name -> Formula -> Unknown */}
                {compound.name && compound.name !== "GeneratedCompound" 
                  ? compound.name 
                  : compound.formula || "Unknown Compound"}
              </h1>
            </div>
          </div>

          <div className="bg-gray-100 dark:bg-slate-700 p-3 rounded-lg mb-6 font-mono text-sm text-gray-600 dark:text-gray-300 break-all border border-gray-200 dark:border-slate-600">
            {compound.smiles}
          </div>

          <div className="grid grid-cols-2 gap-6 mb-8">
            <div className="p-4 bg-gray-50 dark:bg-slate-800 rounded-xl border border-border">
              <p className="text-sm text-muted mb-1">Formula</p>
              <p className="text-lg font-bold text-main">{compound.formula}</p>
            </div>
            <div className="p-4 bg-gray-50 dark:bg-slate-800 rounded-xl border border-border">
              <p className="text-sm text-muted mb-1">Mol Weight</p>
              <p className="text-lg font-bold text-main">{compound.molecular_weight} g/mol</p>
            </div>
            <div className="p-4 bg-gray-50 dark:bg-slate-800 rounded-xl border border-border col-span-2">
              <div className="flex justify-between items-center mb-1">
                <p className="text-sm text-muted">Tg Prediction</p>
                <span className={`text-xs px-2 py-0.5 rounded ${compound.tg > 50 ? 'bg-red-100 text-red-600' : 'bg-green-100 text-green-600'}`}>
                  {compound.tg > 0 ? "High Heat" : "Flexible"}
                </span>
              </div>
              <p className="text-2xl font-bold text-main mb-1">
                {compound.tg} °C
              </p>
              <p className="text-xs text-gray-500 italic">
                "{compound.tg_justification || 'Prediksi berdasarkan struktur kimia.'}"
              </p>
            </div>
          </div>

          <div className="mb-8">
            <h3 className="font-bold text-main mb-2">AI Analysis</h3>
            <p className="text-gray-600 dark:text-gray-300 leading-relaxed text-sm">
              {compound.justifikasi}
            </p>
          </div>

          <button 
            onClick={handleSave}
            disabled={isSaving}
            className="mt-auto w-full bg-slate-900 hover:bg-slate-800 dark:bg-blue-600 dark:hover:bg-blue-700 text-white py-4 rounded-xl font-bold text-lg shadow-lg transform transition hover:-translate-y-1 active:scale-95"
          >
            {isSaving ? "Saving..." : "Save to Library"}
          </button>
        </div>
      </div>

      <h2 className="text-2xl font-bold text-main mt-12 mb-6">Similar Compounds Found</h2>
      
      {/* --- BAGIAN SIMILAR COMPOUNDS (DIPERBAIKI) --- */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        {/* 3. TIDAK PERLU 'any' LAGI KARENA SUDAH ADA INTERFACE */}
        {similar.map((item, idx) => (
          <div key={idx} className="bg-card p-5 rounded-xl border border-border hover:shadow-lg transition-all group">
            <div className="flex justify-between items-center mb-4">
              <span className="bg-green-100 text-green-700 text-xs font-bold px-2 py-1 rounded">
                {item.similarity_percent.toFixed(1)}% Match
              </span>
            </div>
            
            <div className="h-32 bg-white rounded-lg border border-gray-100 mb-4 p-2 flex items-center justify-center">
               <img src={`http://127.0.0.1:8000${item.image_url}`} alt={item.name} className="max-h-full max-w-full object-contain" />
            </div>

            <h3 className="font-bold text-lg text-main mb-1 line-clamp-1" title={item.name}>
              {item.name || "Unknown"}
            </h3>
            <p className="text-xs text-muted mb-4 font-mono truncate">{item.formula}</p>

            <div className="flex justify-between text-sm border-t border-border pt-3">
              <div>
                <p className="text-xs text-muted">MW</p>
                <p className="font-semibold text-main">{item.molecular_weight}</p>
              </div>
              <div className="text-right">
                <p className="text-xs text-muted">Tg</p>
                <p className="font-semibold text-main">{item.tg} °C</p>
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

export default ChemicalDetailPage;