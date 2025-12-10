import { Link } from 'react-router-dom';

function HomePage() {
  return (
    <div className="max-w-5xl mx-auto mt-12 px-4">
      
      {/* Header Section */}
      <div className="text-center mb-12">
        <h1 className="text-4xl md:text-5xl font-extrabold text-gray-900 tracking-tight mb-4">
          Plastic Discovery Agent
        </h1>
        <p className="text-lg text-gray-600 max-w-2xl mx-auto leading-relaxed">
          Masukkan monomer, polimer, atau struktur SMILES untuk memprediksi properti kimia dan jalur sintesis berbasis AI.
        </p>
      </div>

      {/* SEARCH BOX (Desain Baru - Lebih Bersih) */}
      <div className="bg-white p-2 rounded-2xl shadow-lg border border-gray-200 flex items-center max-w-3xl mx-auto mb-20 transition-shadow focus-within:shadow-xl focus-within:border-blue-300"> 
        {/* Input Utama */}
        <input 
          type="text" 
          placeholder="Cari: 'Polyethylene', 'C=C', atau paste SMILES..."
          className="flex-1 px-5 py-3 text-gray-700 placeholder-gray-400 outline-none text-lg bg-transparent"
        />
        
        {/* Tombol Search */}
        <Link 
          to="/explore"
          className="bg-gray-900 hover:bg-blue-600 text-white px-8 py-3 rounded-xl font-bold transition-all duration-300 shadow-md hover:shadow-lg transform hover:-translate-y-0.5"
        >
          SEARCH
        </Link>
      </div>

      {/* MODULE CARDS (Professional & Clean) */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        
        {/* Card 1: Biodegradability */}
        <div className="bg-white p-6 rounded-2xl border border-gray-200 hover:border-blue-500 hover:shadow-lg transition-all cursor-pointer group">
          {/* Icon Container: Abu-abu jadi Biru saat Hover */}
          <div className="w-14 h-14 rounded-xl bg-gray-50 text-gray-600 flex items-center justify-center mb-5 group-hover:bg-blue-600 group-hover:text-white transition-all duration-300">
            {/* Icon Daun (Leaf) */}
            <svg className="w-7 h-7" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M4.318 6.318a4.5 4.5 0 000 6.364L12 20.364l7.682-7.682a4.5 4.5 0 00-6.364-6.364L12 7.636l-1.318-1.318a4.5 4.5 0 00-6.364 0z" /></svg>
          </div>
          <h3 className="text-lg font-bold text-gray-900 mb-2 group-hover:text-blue-600 transition-colors">Biodegradability</h3>
          <p className="text-sm text-gray-500 leading-relaxed">
            Prediksi potensi penguraian alami struktur polimer berdasarkan gugus fungsi.
          </p>
        </div>

        {/* Card 2: Thermal Properties */}
        <div className="bg-white p-6 rounded-2xl border border-gray-200 hover:border-blue-500 hover:shadow-lg transition-all cursor-pointer group">
          <div className="w-14 h-14 rounded-xl bg-gray-50 text-gray-600 flex items-center justify-center mb-5 group-hover:bg-blue-600 group-hover:text-white transition-all duration-300">
            {/* Icon Api/Thermometer */}
            <svg className="w-7 h-7" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M17.657 18.657A8 8 0 016.343 7.343S7 9 9 10c0-2 .5-5 2.986-7C14 5 16.09 5.777 17.656 7.343A7.975 7.975 0 0120 13a7.975 7.975 0 01-2.343 5.657z" /><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9.879 16.121A3 3 0 1012.015 11L11 14H9c0 .768.293 1.536.879 2.121z" /></svg>
          </div>
          <h3 className="text-lg font-bold text-gray-900 mb-2 group-hover:text-blue-600 transition-colors">Thermal Properties</h3>
          <p className="text-sm text-gray-500 leading-relaxed">
            Analisis titik leleh (Tm) dan transisi gelas (Tg) untuk ketahanan panas.
          </p>
        </div>

        {/* Card 3: Synthesis Path */}
        <div className="bg-white p-6 rounded-2xl border border-gray-200 hover:border-blue-500 hover:shadow-lg transition-all cursor-pointer group">
          <div className="w-14 h-14 rounded-xl bg-gray-50 text-gray-600 flex items-center justify-center mb-5 group-hover:bg-blue-600 group-hover:text-white transition-all duration-300">
            {/* Icon Lab Flask */}
            <svg className="w-7 h-7" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" /></svg>
          </div>
          <h3 className="text-lg font-bold text-gray-900 mb-2 group-hover:text-blue-600 transition-colors">Synthesis Path</h3>
          <p className="text-sm text-gray-500 leading-relaxed">
            Rekomendasi jalur reaksi monomer dan katalis untuk sintesis polimer.
          </p>
        </div>

      </div>

    </div>
  );
}

export default HomePage;