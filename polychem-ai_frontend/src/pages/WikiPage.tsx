function WikiPage() {
  return (
    <div className="max-w-6xl mx-auto">
      {/* Header Section */}
      <div className="mb-10 text-center md:text-left">
        <h1 className="text-3xl font-extrabold text-gray-900">Dokumentasi & Panduan</h1>
        <p className="text-gray-500 mt-2">Panduan lengkap untuk memaksimalkan penggunaan PolyChem AI.</p>
      </div>

      {/* Main Grid Layout */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        
        {/* --- CARD 1: CARA PENCARIAN (Full Width) --- */}
        <div className="col-span-1 md:col-span-2 bg-white rounded-2xl border border-gray-200 shadow-sm p-8 hover:shadow-md transition-shadow">
          <div className="flex items-start gap-5">
            {/* Icon Container */}
            <div className="w-12 h-12 rounded-xl bg-blue-50 flex items-center justify-center flex-shrink-0 text-blue-600">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" /></svg>
            </div>
            
            <div className="flex-1">
              <h2 className="text-xl font-bold text-gray-900 mb-3">Cara Memulai Pencarian</h2>
              <p className="text-gray-600 mb-6 leading-relaxed">
                PolyChem AI dirancang fleksibel. Anda tidak perlu menghafal rumus rumit. Cukup gunakan salah satu dari dua metode di bawah ini pada kolom pencarian utama:
              </p>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {/* Metode 1 */}
                <div className="border border-gray-100 bg-gray-50 rounded-lg p-4">
                  <div className="flex items-center gap-2 mb-2">
                    <span className="bg-blue-600 text-white text-xs font-bold px-2 py-0.5 rounded">OPSI 1</span>
                    <span className="font-semibold text-gray-800">Nama Umum</span>
                  </div>
                  <p className="text-sm text-gray-600">Masukkan nama pasar atau nama kimia umum.</p>
                  <div className="mt-2 text-xs font-mono text-blue-600 bg-blue-50 inline-block px-2 py-1 rounded border border-blue-100">
                    Contoh: "Polyethylene", "PET"
                  </div>
                </div>

                {/* Metode 2 */}
                <div className="border border-gray-100 bg-gray-50 rounded-lg p-4">
                  <div className="flex items-center gap-2 mb-2">
                    <span className="bg-blue-600 text-white text-xs font-bold px-2 py-0.5 rounded">OPSI 2</span>
                    <span className="font-semibold text-gray-800">Format SMILES</span>
                  </div>
                  <p className="text-sm text-gray-600">Masukkan string notasi kimia untuk struktur spesifik.</p>
                  <div className="mt-2 text-xs font-mono text-blue-600 bg-blue-50 inline-block px-2 py-1 rounded border border-blue-100">
                    Contoh: CCO, c1ccccc1
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* --- CARD 2: PENJELASAN SMILES (Half Width) --- */}
        <div className="bg-white rounded-2xl border border-gray-200 shadow-sm p-8 hover:shadow-md transition-shadow flex flex-col">
          <div className="flex items-center gap-4 mb-4">
            <div className="w-10 h-10 rounded-lg bg-blue-50 flex items-center justify-center text-blue-600">
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" /></svg>
            </div>
            <h2 className="text-lg font-bold text-gray-900">Apa itu Format SMILES?</h2>
          </div>
          
          <p className="text-gray-600 text-sm mb-4 leading-relaxed flex-grow">
            <strong>SMILES</strong> <em>(Simplified Molecular Input Line Entry System)</em> adalah standar internasional untuk menuliskan struktur 3D molekul menjadi baris teks ASCII.
          </p>

          <div className="bg-slate-100 rounded-xl p-4 text-slate-900 font-mono text-xs border border-slate-100">
            <div className="flex justify-between mb-2 border-b border-slate-700 pb-2">
              <span>Molekul</span>
              <span>SMILES Code</span>
            </div>
            <div className="space-y-1">
              <div className="flex justify-between">
                <span className="text-black">Ethanol</span>
                <span className="text-green-400">CCO</span>
              </div>
              <div className="flex justify-between">
                <span className="text-black">Benzene</span>
                <span className="text-green-400">c1ccccc1</span>
              </div>
              <div className="flex justify-between">
                <span className="text-black">Carbon Dioxide</span>
                <span className="text-green-400">O=C=O</span>
              </div>
            </div>
          </div>
        </div>

        {/* --- CARD 3: INTERPRETASI SKOR (Half Width) --- */}
        <div className="bg-white rounded-2xl border border-gray-200 shadow-sm p-8 hover:shadow-md transition-shadow flex flex-col">
          <div className="flex items-center gap-4 mb-4">
            <div className="w-10 h-10 rounded-lg bg-green-50 flex items-center justify-center text-green-600">
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" /></svg>
            </div>
            <h2 className="text-lg font-bold text-gray-900">AI Confidence Score</h2>
          </div>

          <p className="text-gray-600 text-sm mb-6 leading-relaxed">
            Seberapa yakin AI terhadap prediksinya? Gunakan panduan ini untuk memverifikasi hasil.
          </p>

          <div className="space-y-3">
            {/* Level High */}
            <div className="flex items-center gap-3 p-3 rounded-lg border border-green-100 bg-green-50/50">
              <div className="w-2 h-10 rounded-full bg-green-500"></div>
              <div>
                <p className="text-xs font-bold text-green-700 uppercase tracking-wider">90% - 100%</p>
                <p className="text-xs text-gray-600 font-medium">Sangat Akurat (Data validasi tersedia)</p>
              </div>
            </div>

            {/* Level Medium */}
            <div className="flex items-center gap-3 p-3 rounded-lg border border-yellow-100 bg-yellow-50/50">
              <div className="w-2 h-10 rounded-full bg-yellow-500"></div>
              <div>
                <p className="text-xs font-bold text-yellow-700 uppercase tracking-wider">70% - 89%</p>
                <p className="text-xs text-gray-600 font-medium">Moderat (Disarankan verifikasi lab)</p>
              </div>
            </div>

            {/* Level Low */}
            <div className="flex items-center gap-3 p-3 rounded-lg border border-red-100 bg-red-50/50">
              <div className="w-2 h-10 rounded-full bg-red-500"></div>
              <div>
                <p className="text-xs font-bold text-red-700 uppercase tracking-wider">&lt; 70%</p>
                <p className="text-xs text-gray-600 font-medium">Rendah (Struktur kimia tidak umum)</p>
              </div>
            </div>
          </div>
        </div>

      </div>
    </div>
  );
}

export default WikiPage;