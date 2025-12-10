import { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import { auth } from '../lib/firebase';
import { getLibrary, removeFromLibrary } from '../services/dbService'; 
import { onAuthStateChanged } from 'firebase/auth';
import { type ChemicalData } from '../services/mockData';

// Gabungkan tipe data Chemical dengan savedAt
type LibraryItem = ChemicalData & { savedAt: string };

function LibraryPage() {
  const [libraryData, setLibraryData] = useState<LibraryItem[]>([]);
  const [loading, setLoading] = useState(true);
  const [userUid, setUserUid] = useState<string | null>(null);

  // --- 1. DEFINISIKAN FUNGSI FETCH DATA DULU (SUPAYA BISA DIPANGGIL) ---
  const fetchData = async (uid: string) => {
    const data = await getLibrary(uid);
    setLibraryData(data as LibraryItem[]);
    setLoading(false);
  };

  // --- 2. BARU PANGGIL DI USE EFFECT ---
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, async (user) => {
      if (user) {
        setUserUid(user.uid);
        fetchData(user.uid); // Sekarang aman, karena fetchData sudah ada di atas
      } else {
        setLibraryData([]);
        setLoading(false);
      }
    });
    return () => unsubscribe();
  }, []);

  // --- FUNGSI HAPUS ---
  const handleDelete = async (e: React.MouseEvent, id: number) => {
    e.preventDefault(); 
    if (!userUid) return;

    if (window.confirm("Yakin ingin menghapus item ini dari library?")) {
      // Optimistic Update
      setLibraryData((prev) => prev.filter((item) => item.id !== id));

      const success = await removeFromLibrary(userUid, id);
      if (!success) {
        alert("Gagal menghapus data.");
        fetchData(userUid); // Fetch ulang jika gagal
      }
    }
  };

  return (
    <div className="max-w-7xl mx-auto">
      <div className="mb-8">
        <h1 className="text-2xl font-bold text-gray-900">My Library</h1>
        <p className="text-gray-600">Koleksi senyawa yang tersimpan.</p>
      </div>

      {loading ? (
        <div className="text-center py-10 text-gray-500">Memuat library...</div>
      ) : libraryData.length > 0 ? (
        <div className="grid grid-cols-1 gap-6 sm:grid-cols-2 lg:grid-cols-3">
          {libraryData.map((item, index) => (
            <Link 
              to={`/chemical/${item.id}`} 
              key={index} 
              className="bg-white block overflow-hidden shadow-sm rounded-xl border border-gray-200 hover:border-blue-400 hover:shadow-md transition-all relative group"
            >
              <div className="p-6">
                <div className="flex justify-between items-start mb-4">
                  <div className="inline-flex items-center justify-center w-10 h-10 rounded-lg bg-blue-50 text-blue-700 font-bold text-lg">
                    {item.name ? item.name.charAt(0) : '?'}
                  </div>
                  
                  {/* TOMBOL DELETE */}
                  <button 
                    onClick={(e) => handleDelete(e, item.id)}
                    className="p-2 text-gray-400 hover:text-red-600 hover:bg-red-50 rounded-full transition-colors z-10"
                    title="Hapus"
                  >
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" /></svg>
                  </button>
                </div>
                
                <h3 className="text-lg font-bold text-gray-900 mb-1">{item.name}</h3>
                <p className="text-xs text-gray-500 font-mono bg-gray-50 p-2 rounded truncate mb-2">
                  {item.smiles}
                </p>
                <div className="text-xs text-gray-400">
                  Disimpan: {new Date(item.savedAt).toLocaleDateString()}
                </div>
              </div>
            </Link>
          ))}
        </div>
      ) : (
        <div className="p-12 text-center bg-white rounded-xl border border-dashed border-gray-300">
          <p className="text-gray-500 mb-4">Belum ada data di library.</p>
          <Link to="/explore" className="text-blue-600 hover:underline">Cari Senyawa</Link>
        </div>
      )}
    </div>
  );
}

export default LibraryPage;