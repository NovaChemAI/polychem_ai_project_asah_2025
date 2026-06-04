import { useLocation, useNavigate } from "react-router-dom";
import { useEffect, useState } from "react";
import toast from "react-hot-toast";
import { saveToLibrary } from "../services/dbService";
import { auth } from "../lib/firebase";
import { buildAssetUrl, normalizePrediction } from "../services/apiService";

function ChemicalDetailPage() {
  const location = useLocation();
  const navigate = useNavigate();
  const [isSaving, setIsSaving] = useState(false);

  const predictionData = normalizePrediction(location.state?.predictionData);

  useEffect(() => {
    if (!predictionData) {
      navigate("/");
    }
  }, [predictionData, navigate]);

  if (!predictionData) return null;

  const compound = predictionData.new_compound;
  const similar = predictionData.similar_compounds;

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
      category: compound.polymer_class || "Novel Compound",
    };

    const success = await saveToLibrary(user.uid, payload);
    if (success) toast.success("Berhasil disimpan ke Library!");
    else toast.error("Gagal menyimpan.");

    setIsSaving(false);
  };

  return (
    <div className="max-w-6xl mx-auto px-4 py-8">
      <button
        onClick={() => navigate(-1)}
        className="mb-6 text-gray-500 hover:text-blue-600 flex items-center gap-2 transition-colors font-medium"
      >
        ← Kembali
      </button>

      {/* --- KARTU UTAMA (HASIL PREDIKSI) --- */}
      <div className="bg-white dark:bg-slate-900 rounded-3xl shadow-2xl border border-border overflow-hidden flex flex-col md:flex-row transition-all">
        {/* Sisi Kiri: Gambar Struktur Utama */}
        <div className="md:w-1/3 bg-gray-50 dark:bg-slate-800/50 p-10 flex items-center justify-center border-r border-border relative">
          {compound.image_url ? (
            <img
              src={buildAssetUrl(compound.image_url)}
              alt="Structure"
              className="w-full h-auto max-w-[280px] object-contain mix-blend-multiply dark:mix-blend-normal filter dark:invert"
              onError={(e) => {
                (e.target as HTMLImageElement).src =
                  "https://via.placeholder.com/300?text=No+Structure";
              }}
            />
          ) : (
            <div className="text-gray-400 font-medium">No Structure Image</div>
          )}
        </div>

        {/* Sisi Kanan: Detail Informasi */}
        <div className="p-8 md:p-10 md:w-2/3 flex flex-col">
          <div className="mb-6">
            <span className="inline-block px-4 py-1.5 text-[10px] font-black tracking-widest text-blue-700 uppercase bg-blue-100 dark:bg-blue-900/30 dark:text-blue-400 rounded-lg mb-3">
              {compound.polymer_class || "Novel Compound"}
            </span>
            {/* PERBAIKAN NAMA: Sekarang konsisten "Unknown Compound" jika tidak ada nama spesifik */}
            <h1 className="text-3xl md:text-4xl font-black text-slate-900 dark:text-white leading-tight">
              {compound.name && compound.name !== "GeneratedCompound"
                ? compound.name
                : "Unknown Compound"}
            </h1>
          </div>

          {/* SMILES Box */}
          <div className="bg-slate-50 dark:bg-slate-800/80 p-4 rounded-xl mb-8 font-mono text-xs text-slate-500 dark:text-slate-400 break-all border border-slate-200 dark:border-slate-700 relative group">
            <span className="absolute -top-2 left-3 bg-white dark:bg-slate-900 px-2 text-[9px] font-bold text-slate-400 border border-slate-200 dark:border-slate-700 rounded uppercase">
              SMILES Notation
            </span>
            {compound.smiles}
          </div>

          {/* Grid Informasi Utama */}
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 mb-8">
            <div className="p-4 bg-white dark:bg-slate-800 rounded-2xl border border-slate-100 dark:border-slate-700 shadow-sm">
              <p className="text-[10px] font-bold text-slate-400 uppercase tracking-wider mb-1">
                Formula predict
              </p>
              <p className="text-xl font-bold text-slate-800 dark:text-slate-200">
                {compound.formula}
              </p>
            </div>
            <div className="p-4 bg-white dark:bg-slate-800 rounded-2xl border border-slate-100 dark:border-slate-700 shadow-sm">
              <p className="text-[10px] font-bold text-slate-400 uppercase tracking-wider mb-1">
                Mol. Weight predict
              </p>
              <p className="text-xl font-bold text-slate-800 dark:text-slate-200">
                {compound.molecular_weight.toFixed(2)}{" "}
                <span className="text-xs font-normal opacity-60">g/mol</span>
              </p>
            </div>
            <div className="p-5 bg-blue-50/50 dark:bg-blue-900/10 rounded-2xl border border-blue-100 dark:border-blue-800/50 col-span-1 sm:col-span-2">
              <div className="flex justify-between items-center mb-2">
                <p className="text-[10px] font-bold text-blue-500 uppercase tracking-wider">
                  Glass Transition Temperature (Tg) predict
                </p>
                <span
                  className={`text-[10px] px-2 py-1 rounded-md font-black uppercase ${compound.tg > 50 ? "bg-orange-100 text-orange-600" : "bg-green-100 text-green-600"}`}
                >
                  {compound.tg > 50 ? "Heat Resistant" : "Highly Flexible"}
                </span>
              </div>
              <p className="text-4xl font-black text-slate-900 dark:text-white mb-2">
                {compound.tg} °C
              </p>
              <p className="text-xs text-slate-500 dark:text-slate-400 italic bg-white/50 dark:bg-slate-900/50 p-2 rounded-lg border border-blue-50 dark:border-blue-900/20">
                "{compound.tg_justification}"
              </p>
            </div>
          </div>

          {/* AI Analysis */}
          <div className="mb-10">
            <h3 className="text-sm font-black text-slate-900 dark:text-white uppercase tracking-widest mb-3 flex items-center gap-2">
              <span className="w-2 h-2 bg-blue-500 rounded-full"></span>{" "}
              Justification
            </h3>
            <p className="text-slate-600 dark:text-slate-300 leading-relaxed text-sm bg-slate-50 dark:bg-slate-800/30 p-4 rounded-xl border border-slate-100 dark:border-slate-700">
              {compound.justifikasi}
            </p>
          </div>

          <button
            onClick={handleSave}
            disabled={isSaving}
            className="mt-auto w-full bg-slate-900 hover:bg-black dark:bg-blue-600 dark:hover:bg-blue-500 text-white py-4 rounded-2xl font-black text-lg shadow-xl transform transition hover:-translate-y-1 active:scale-95 disabled:opacity-50 disabled:cursor-not-allowed"
          >
            {isSaving ? "PROCESSING..." : "SAVE TO LIBRARY"}
          </button>
        </div>
      </div>

      {/* --- BAGIAN SIMILAR COMPOUNDS --- */}
      <div className="mt-16">
        <div className="flex items-center gap-4 mb-8">
          <h2 className="text-2xl font-black text-slate-900 dark:text-white">
            Similar Compounds Found
          </h2>
          <div className="h-px flex-1 bg-slate-200 dark:bg-slate-800"></div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          {similar.map((item, idx) => (
            <div
              key={idx}
              className="bg-white dark:bg-slate-900 flex flex-col p-6 rounded-2xl border border-slate-200 dark:border-slate-800 hover:shadow-2xl hover:border-blue-300 dark:hover:border-blue-800 transition-all duration-300 group"
            >
              {/* Badge Match */}
              <div className="flex justify-between items-center mb-5">
                <span className="bg-green-100 dark:bg-green-900/30 text-green-700 dark:text-green-400 text-[10px] font-black px-2.5 py-1 rounded-lg uppercase tracking-wider">
                  {item.similarity_percent.toFixed(1)}% Match
                </span>
                <span className="text-slate-300 dark:text-slate-700 font-black text-xs">
                  #{idx + 1}
                </span>
              </div>

              {/* Struktur Kimia */}
              <div className="h-44 bg-slate-50 dark:bg-slate-800/50 rounded-xl border border-slate-100 dark:border-slate-800 mb-6 p-4 flex items-center justify-center overflow-hidden relative">
                <img
                  src={buildAssetUrl(item.image_url)}
                  alt={item.name}
                  className="max-h-full max-w-full object-contain mix-blend-multiply dark:mix-blend-normal filter dark:invert group-hover:scale-110 transition-transform duration-500"
                  onError={(e) => {
                    (e.target as HTMLImageElement).src =
                      "https://via.placeholder.com/150?text=No+Image";
                  }}
                />
              </div>

              {/* Info Senyawa */}
              <div className="mb-5">
                <h3
                  className="font-black text-lg text-slate-900 dark:text-white mb-1 truncate"
                  title={item.name}
                >
                  {item.name || "Unknown Compound"}
                </h3>
                <p className="text-[11px] font-mono text-slate-400 bg-slate-50 dark:bg-slate-800/80 p-2 rounded-lg border border-slate-100 dark:border-slate-700 truncate">
                  {item.formula}
                </p>
              </div>

              {/* Data Grid */}
              <div className="grid grid-cols-2 gap-3 mb-6">
                <div className="p-3 bg-slate-50 dark:bg-slate-800/30 rounded-xl border border-slate-100 dark:border-slate-800">
                  <p className="text-[9px] uppercase font-black text-slate-400 tracking-widest mb-1">
                    Mol. Wt
                  </p>
                  <p className="font-bold text-slate-800 dark:text-slate-200 text-sm">
                    {item.molecular_weight.toFixed(2)}
                  </p>
                </div>
                <div className="p-3 bg-slate-50 dark:bg-slate-800/30 rounded-xl border border-slate-100 dark:border-slate-800">
                  <p className="text-[9px] uppercase font-black text-slate-400 tracking-widest mb-1">
                    Tg Predict
                  </p>
                  <p className="font-bold text-slate-800 dark:text-slate-200 text-sm">
                    {item.tg} °C
                  </p>
                </div>
              </div>

              {/* Analisis Singkat Otomatis */}
              <div className="mt-auto pt-5 border-t border-dashed border-slate-200 dark:border-slate-800">
                <div className="flex items-start gap-2">
                  <div
                    className={`w-1.5 h-1.5 rounded-full mt-1.5 shrink-0 ${item.tg > 50 ? "bg-orange-400" : "bg-blue-400"}`}
                  ></div>
                  <p className="text-[11px] text-slate-500 dark:text-slate-400 leading-relaxed italic">
                    Struktur ini teridentifikasi memiliki kemiripan fungsional
                    tinggi dengan target novel compound, menunjukkan
                    karakteristik termal stabil pada {item.tg} °C.
                  </p>
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

export default ChemicalDetailPage;
