import { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import { chemicalDatabase, type ChemicalData } from '../services/mockData';

// Import Firebase & Service
import { auth } from '../lib/firebase';
import { saveSearchToHistory } from '../services/dbService';
import { onAuthStateChanged } from 'firebase/auth';

function ExplorePage() {
  const [searchTerm, setSearchTerm] = useState('');
  const [userUid, setUserUid] = useState<string | null>(null);

  // 1. Cek Status Login User (Siapa yang sedang mencari?)
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

  // 2. Fungsi Menangani Tombol Enter
  const handleKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter' && searchTerm.trim() !== '') {
      console.log("User menekan Enter:", searchTerm);
      
      // Hanya simpan jika user sudah login
      if (userUid) {
        saveSearchToHistory(userUid, searchTerm);
      }
    }
  };

  // 3. Logika Filter Data (Tetap sama)
  const filteredData = chemicalDatabase.filter((item: ChemicalData) => {
    return item.name.toLowerCase().includes(searchTerm.toLowerCase()) || 
           item.smiles.toLowerCase().includes(searchTerm.toLowerCase());
  });

  return (
    <div className="max-w-7xl mx-auto">
      
      {/* Header & Search Section */}
      <div className="mb-8 flex flex-col md:flex-row md:items-end justify-between gap-4">
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Polymer Analysis</h1>
          <p className="text-gray-600">Database material plastik dan polimer.</p>
        </div>

        {/* Search Input */}
        <div className="relative w-full md:w-96">
          <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
            <svg className="h-5 w-5 text-gray-400" viewBox="0 0 20 20" fill="currentColor">
              <path fillRule="evenodd" d="M8 4a4 4 0 100 8 4 4 0 000-8zM2 8a6 6 0 1110.89 3.476l4.817 4.817a1 1 0 01-1.414 1.414l-4.816-4.816A6 6 0 012 8z" clipRule="evenodd" />
            </svg>
          </div>
          <input
            type="text"
            className="block w-full pl-10 pr-3 py-2 border border-gray-300 rounded-lg leading-5 bg-white placeholder-gray-500 focus:outline-none focus:placeholder-gray-400 focus:border-blue-300 focus:ring focus:ring-blue-200 sm:text-sm transition duration-150 ease-in-out"
            placeholder="Cari polimer (cth: PET, PLA)..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            onKeyDown={handleKeyDown} // <--- Event Handler Baru
          />
        </div>
      </div>

      {/* Grid Hasil */}
      <div className="grid grid-cols-1 gap-6 sm:grid-cols-2 lg:grid-cols-3">
        
        {filteredData.length > 0 ? (
          filteredData.map((item) => (
            <Link 
              to={`/chemical/${item.id}`} 
              key={item.id} 
              className="bg-white block overflow-hidden shadow-sm rounded-xl border border-gray-200 hover:border-blue-400 hover:shadow-md transition-all"
            >
              <div className="p-6">
                <div className="flex justify-between items-start mb-4">
                  <div className="inline-flex items-center justify-center w-10 h-10 rounded-lg bg-gray-100 text-gray-700 font-bold text-lg">
                    {item.name.charAt(0)}
                  </div>
                  <span className={`px-3 py-1 rounded-full text-xs font-medium ${
                    item.category === 'Biodegradable' ? 'bg-green-100 text-green-800' : 
                    item.category === 'Thermoplastic' ? 'bg-blue-100 text-blue-800' : 
                    item.category === 'Thermoset' ? 'bg-purple-100 text-purple-800' :
                    'bg-gray-100 text-gray-800'
                  }`}>
                    {item.category}
                  </span>
                </div>
                
                <h3 className="text-lg font-bold text-gray-900 mb-1">{item.name}</h3>
                
                <div className="mb-4">
                   <p className="text-xs text-gray-400 mb-1">Structure (SMILES):</p>
                   <p className="text-xs text-gray-600 font-mono bg-gray-50 p-2 rounded border border-gray-100 break-all">
                    {item.smiles}
                   </p>
                </div>
                
                <div className="border-t border-gray-100 pt-3">
                  <p className="text-xs text-gray-400">Key Properties:</p>
                  <p className="text-sm text-gray-700 font-medium mt-1 line-clamp-2">
                    {item.properties}
                  </p>
                </div>
                
                <div className="flex items-center justify-between pt-4 mt-2">
                  <span className="text-xs text-gray-500">Prediction Score:</span>
                  <span className="text-sm font-bold text-blue-600">{item.score}%</span>
                </div>
              </div>
            </Link>
          ))
        ) : (
          <div className="col-span-full text-center py-12 bg-white rounded-xl border border-dashed border-gray-300">
            <p className="text-gray-500">Tidak ada data polimer yang cocok untuk "{searchTerm}".</p>
          </div>
        )}
        
      </div>
    </div>
  );
}

export default ExplorePage;