import { useState, useEffect } from 'react';
import { useParams, Link } from 'react-router-dom';
import { chemicalDatabase } from '../services/mockData';

// Import Logic
import { auth } from '../lib/firebase';
import { saveToLibrary } from '../services/dbService';
import { onAuthStateChanged } from 'firebase/auth';

function ChemicalDetailPage() {
  const { id } = useParams();
  const chemical = chemicalDatabase.find(item => item.id === Number(id));
  
  const [userUid, setUserUid] = useState<string | null>(null);
  const [isSaved, setIsSaved] = useState(false);
  const [isSaving, setIsSaving] = useState(false);

  // 1. Cek Login
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (user) => {
      if (user) {
        setUserUid(user.uid);
      } else {
        setUserUid(null);
      }
    });
    return () => unsubscribe();
  }, []);

  // 2. Handle Save
  const handleSave = async () => {
    // Validasi: Harus Login
    if (!userUid) {
      alert("Silakan login terlebih dahulu untuk menyimpan ke Library.");
      return;
    }
    
    if (chemical) {
      setIsSaving(true);
      // Panggil fungsi service yang sudah diperbaiki
      const success = await saveToLibrary(userUid, chemical);
      
      if (success) {
        setIsSaved(true);
      } else {
        alert("Gagal menyimpan. Cek koneksi internet.");
      }
      setIsSaving(false);
    }
  };

  if (!chemical) {
    return (
      <div className="p-20 text-center flex flex-col items-center justify-center">
        <div className="text-6xl mb-4">🧪</div>
        <h2 className="text-2xl font-bold text-gray-900">Senyawa Tidak Ditemukan</h2>
        <Link to="/explore" className="bg-gray-900 text-white px-6 py-2 mt-4 rounded-lg">Kembali</Link>
      </div>
    );
  }

  return (
    <div className="max-w-5xl mx-auto">
      <div className="mb-6 text-sm text-gray-500 flex items-center gap-2">
        <Link to="/explore" className="hover:text-blue-600 hover:underline">Explore</Link>
        <span>/</span>
        <span className="font-semibold text-gray-900">{chemical.name}</span>
      </div>

      <div className="bg-white rounded-2xl shadow-sm border border-gray-200 overflow-hidden">
        <div className="p-8 md:flex gap-10">
          
          {/* Gambar Kiri */}
          <div className="w-full md:w-1/3">
            <div className="aspect-square bg-white rounded-2xl border border-slate-200 flex flex-col items-center justify-center p-4 relative group overflow-hidden">
               {chemical.image ? (
                 <img src={chemical.image} alt={chemical.name} className="w-full h-full object-contain p-2" />
               ) : (
                 <div className="text-8xl z-10 mb-4">⚗️</div>
               )}
            </div>
          </div>

          {/* Info Kanan */}
          <div className="w-full md:w-2/3 mt-6 md:mt-0 flex flex-col justify-center">
            <div className="flex justify-between items-start">
              <div>
                <h1 className="text-3xl font-extrabold text-gray-900 mb-2">{chemical.name}</h1>
                <span className="inline-block px-3 py-1 rounded-full text-sm font-medium bg-blue-100 text-blue-800">
                  {chemical.category}
                </span>
              </div>
              <div className="text-right bg-blue-50 px-4 py-2 rounded-lg border border-blue-100">
                <p className="text-xs text-blue-600 font-semibold">AI Confidence</p>
                <p className="text-3xl font-black text-blue-700">{chemical.score}%</p>
              </div>
            </div>

            <div className="mt-8 space-y-6">
               <div>
                 <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wide mb-1">SMILES</h3>
                 <div className="font-mono text-sm bg-slate-100 text-slate-600 p-2 rounded border border-slate-200 break-all">
                   {chemical.smiles}
                 </div>
               </div>
               <div>
                 <h3 className="text-xs font-bold text-gray-400 uppercase tracking-wide mb-2">Properties</h3>
                 <div className="bg-gray-50 p-4 rounded-xl border border-gray-100 text-gray-700 text-sm">
                    {chemical.properties}
                 </div>
               </div>

               <div className="flex gap-3 pt-4">
                  {/* TOMBOL SAVE */}
                  <button 
                    onClick={handleSave}
                    disabled={isSaved || isSaving}
                    className={`px-6 py-3 rounded-xl text-sm font-semibold transition-colors border flex items-center gap-2 ${
                      isSaved 
                        ? 'bg-green-50 text-green-700 border-green-200 cursor-default' 
                        : 'bg-white text-gray-700 border-gray-300 hover:bg-gray-50'
                    }`}
                  >
                    {isSaving ? 'Saving...' : isSaved ? 'Saved to Library' : 'Save to Project'}
                  </button>
               </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

export default ChemicalDetailPage;