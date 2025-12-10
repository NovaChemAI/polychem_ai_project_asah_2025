import { useState, useEffect } from 'react';
import { useParams, Link } from 'react-router-dom';
import { chemicalDatabase } from '../services/mockData';

// Import Logic
import { auth } from '../lib/firebase';
import { saveToLibrary, removeFromLibrary, checkIsSaved } from '../services/dbService'; // Tambah fungsi baru
import { onAuthStateChanged } from 'firebase/auth';

function ChemicalDetailPage() {
  const { id } = useParams();
  const chemical = chemicalDatabase.find(item => item.id === Number(id));
  
  const [userUid, setUserUid] = useState<string | null>(null);
  const [isSaved, setIsSaved] = useState(false); // Status simpan
  const [loading, setLoading] = useState(false); // Loading tombol

  // 1. Cek Login & Cek Apakah Sudah Disimpan
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, async (user) => {
      if (user) {
        setUserUid(user.uid);
        if (chemical) {
          // Cek database apakah item ini sudah ada
          const saved = await checkIsSaved(user.uid, chemical.id);
          setIsSaved(saved);
        }
      } else {
        setUserUid(null);
      }
    });
    return () => unsubscribe();
  }, [chemical]); 

  // 2. Handle Toggle (Save / Unsave)
  const handleToggleSave = async () => {
    if (!userUid) {
      alert("Silakan login terlebih dahulu.");
      return;
    }
    
    if (!chemical) return;

    setLoading(true);

    if (isSaved) {
      // --- LOGIKA HAPUS ---
      const success = await removeFromLibrary(userUid, chemical.id);
      if (success) {
        setIsSaved(false); // Ubah status jadi belum disimpan
      } else {
        alert("Gagal menghapus.");
      }
    } else {
      // --- LOGIKA SIMPAN ---
      const success = await saveToLibrary(userUid, chemical);
      if (success) {
        setIsSaved(true); // Ubah status jadi tersimpan
      } else {
        alert("Gagal menyimpan.");
      }
    }
    setLoading(false);
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
                 {/* TOMBOL TOGGLE (SAVE / UNSAVE) */}
                 <button 
                   onClick={handleToggleSave}
                   disabled={loading}
                   className={`px-6 py-3 rounded-xl text-sm font-semibold transition-colors border flex items-center gap-2 ${
                     isSaved 
                       ? 'bg-red-50 text-red-600 border-red-200 hover:bg-red-100' // Tampilan UNSAVE
                       : 'bg-gray-900 text-white border-transparent hover:bg-gray-800' // Tampilan SAVE
                   } ${loading ? 'opacity-50 cursor-not-allowed' : ''}`}
                 >
                   {loading ? 'Processing...' : isSaved ? (
                     <>
                       {/* Ikon Sampah */}
                       <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" /></svg>
                       Remove from Library
                     </>
                   ) : (
                     <>
                       {/* Ikon Bookmark */}
                       <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" /></svg>
                       Save to Library
                     </>
                   )}
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