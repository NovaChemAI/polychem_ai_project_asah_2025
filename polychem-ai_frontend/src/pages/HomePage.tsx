// src/pages/HomePage.tsx
import { Link } from 'react-router-dom';

function HomePage() {
  return (
    <div className="max-w-5xl mx-auto mt-10">
      
      {/* Header Section */}
      <div className="text-center mb-12">
        <h1 className="text-4xl font-extrabold text-gray-900 tracking-tight mb-4">
          Plastic Discovery Agent
        </h1>
        <p className="text-lg text-gray-600">
          Masukkan monomer, polimer, atau struktur SMILES untuk prediksi properti berbasis AI.
        </p>
      </div>

      {/* ASKCOS-style Search Box */}
      <div className="bg-white p-2 rounded-2xl shadow-xl border border-gray-200 flex items-center max-w-3xl mx-auto mb-16">
        {/* Tombol Draw (Fitur masa depan) */}
        <button className="p-3 text-gray-400 hover:text-gray-600 border-r border-gray-200 mr-2" title="Draw Structure">
           ✏️ <span className="text-xs font-bold ml-1">DRAW</span>
        </button>
        
        {/* Input Utama */}
        <input 
          type="text" 
          placeholder="Cari: 'Polyethylene', 'C=C', atau paste SMILES..."
          className="flex-1 px-4 py-3 text-gray-700 placeholder-gray-400 outline-none text-lg"
        />
        
        {/* Tombol Search */}
        <Link 
          to="/explore"
          className="bg-gray-900 hover:bg-gray-800 text-white px-8 py-3 rounded-xl font-bold transition-colors"
        >
          SEARCH
        </Link>
      </div>

      {/* Module Cards (Quick Access) */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        
        {/* Card 1 */}
        <div className="bg-white p-6 rounded-xl border border-gray-200 hover:border-blue-500 hover:shadow-md transition-all cursor-pointer group">
          <div className="w-12 h-12 bg-blue-50 rounded-lg flex items-center justify-center mb-4 group-hover:bg-blue-100 transition-colors">
            <span className="text-2xl">♻️</span>
          </div>
          <h3 className="font-bold text-gray-900 mb-2">Biodegradability</h3>
          <p className="text-sm text-gray-500">Prediksi apakah struktur polimer dapat terurai secara alami.</p>
        </div>

        {/* Card 2 */}
        <div className="bg-white p-6 rounded-xl border border-gray-200 hover:border-purple-500 hover:shadow-md transition-all cursor-pointer group">
          <div className="w-12 h-12 bg-purple-50 rounded-lg flex items-center justify-center mb-4 group-hover:bg-purple-100 transition-colors">
            <span className="text-2xl">🔥</span>
          </div>
          <h3 className="font-bold text-gray-900 mb-2">Thermal Properties</h3>
          <p className="text-sm text-gray-500">Analisis titik leleh (Tm) dan transisi gelas (Tg) material.</p>
        </div>

        {/* Card 3 */}
        <div className="bg-white p-6 rounded-xl border border-gray-200 hover:border-green-500 hover:shadow-md transition-all cursor-pointer group">
           <div className="w-12 h-12 bg-green-50 rounded-lg flex items-center justify-center mb-4 group-hover:bg-green-100 transition-colors">
            <span className="text-2xl">⚗️</span>
          </div>
          <h3 className="font-bold text-gray-900 mb-2">Synthesis Path</h3>
          <p className="text-sm text-gray-500">Temukan jalur reaksi monomer untuk membentuk polimer target.</p>
        </div>

      </div>

    </div>
  );
}

export default HomePage;