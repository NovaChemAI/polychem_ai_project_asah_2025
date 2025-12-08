import { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import { auth } from '../lib/firebase';
import { getLibrary } from '../services/dbService';
import { onAuthStateChanged } from 'firebase/auth';
import { type ChemicalData } from '../services/mockData';

// Gabungkan tipe data Chemical dengan savedAt
type LibraryItem = ChemicalData & { savedAt: string };

function LibraryPage() {
  const [libraryData, setLibraryData] = useState<LibraryItem[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, async (user) => {
      if (user) {
        const data = await getLibrary(user.uid);
        setLibraryData(data as LibraryItem[]);
      } else {
        // Jika user tidak login, kosongkan data
        setLibraryData([]);
      }
      setLoading(false);
    });
    return () => unsubscribe();
  }, []);

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
              className="bg-white block overflow-hidden shadow-sm rounded-xl border border-gray-200 hover:border-blue-400 hover:shadow-md transition-all"
            >
              <div className="p-6">
                <div className="flex justify-between items-start mb-4">
                  <div className="inline-flex items-center justify-center w-10 h-10 rounded-lg bg-blue-50 text-blue-700 font-bold text-lg">
                    {item.name ? item.name.charAt(0) : '?'}
                  </div>
                  <span className="bg-gray-100 text-gray-600 px-2 py-1 rounded text-xs">Saved</span>
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