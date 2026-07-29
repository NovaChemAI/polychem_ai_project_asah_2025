function WikiPage() {
  return (
    <div className="max-w-6xl mx-auto pb-10">
      {/* Header Section */}
      <div className="mb-10 text-center md:text-left pt-6">
        <h1 className="text-3xl font-extrabold text-main">Dokumentasi & Panduan</h1>
        <p className="text-muted mt-2">Panduan lengkap untuk memaksimalkan penggunaan PolyChem AI.</p>
      </div>

      {/* Main Grid Layout */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">

        {/* --- CARD 1: CARA PENCARIAN (Full Width) --- */}
        <div className="col-span-1 md:col-span-2 bg-card rounded-2xl border border-border shadow-sm p-8 hover:shadow-md transition-shadow">
          <div className="flex items-start gap-5">
            <div className="w-12 h-12 rounded-xl bg-blue-50 dark:bg-blue-900/30 flex items-center justify-center flex-shrink-0 text-blue-600 dark:text-blue-400">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" /></svg>
            </div>

            <div className="flex-1">
              <h2 className="text-xl font-bold text-main mb-3">Cara Memulai Pencarian</h2>
              <p className="text-muted mb-6 leading-relaxed">
                Saat ini, PolyChem secara spesifik dirancang untuk menganalisis struktur kimia menggunakan format standar industri. Untuk mendapatkan hasil prediksi yang akurat, mohon gunakan format input berikut:
              </p>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div className="col-span-2 border border-blue-100 dark:border-blue-900 bg-blue-50/50 dark:bg-blue-900/10 rounded-lg p-6 flex flex-col md:flex-row items-center gap-6">
                  <div className="flex-1">
                    <div className="flex items-center gap-2 mb-2">
                       <span className="bg-blue-600 text-white text-xs font-bold px-2 py-0.5 rounded">REQUIRED</span>
                       <span className="font-semibold text-main">Format SMILES</span>
                    </div>
                    <p className="text-sm text-muted mb-3">
                        Masukkan string notasi kimia (SMILES) yang merepresentasikan struktur molekul polimer atau monomer Anda.
                    </p>
                    <div className="flex flex-wrap gap-2">
                        <div className="text-xs font-mono text-blue-700 dark:text-blue-300 bg-white dark:bg-slate-800 px-3 py-1.5 rounded border border-blue-200 dark:border-slate-600 shadow-sm">
                            CCO (Ethanol)
                        </div>
                        <div className="text-xs font-mono text-blue-700 dark:text-blue-300 bg-white dark:bg-slate-800 px-3 py-1.5 rounded border border-blue-200 dark:border-slate-600 shadow-sm">
                            C=C (Ethylene)
                        </div>
                        <div className="text-xs font-mono text-blue-700 dark:text-blue-300 bg-white dark:bg-slate-800 px-3 py-1.5 rounded border border-blue-200 dark:border-slate-600 shadow-sm">
                            c1ccccc1 (Benzene)
                        </div>
                    </div>
                  </div>

                  <div className="hidden md:block w-px h-24 bg-blue-200 dark:bg-blue-800 mx-4"></div>

                  <div className="md:w-1/3 text-xs text-gray-500 dark:text-gray-400 italic text-center md:text-left">
                    "Pastikan penulisan huruf besar dan kecil sesuai (Case Sensitive) untuk atom seperti Cl (Chlorine) atau Br (Bromine)."
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* --- CARD 2: PENJELASAN SMILES (Half Width) --- */}
        <div className="bg-card rounded-2xl border border-border shadow-sm p-8 hover:shadow-md transition-shadow flex flex-col">
          <div className="flex items-center gap-4 mb-4">
            <div className="w-10 h-10 rounded-lg bg-blue-50 dark:bg-blue-900/30 flex items-center justify-center text-blue-600 dark:text-blue-400">
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" /></svg>
            </div>
            <h2 className="text-lg font-bold text-main">Apa itu Format SMILES?</h2>
          </div>

          <p className="text-muted text-sm mb-4 leading-relaxed">
            <strong>SMILES</strong> <em>(Simplified Molecular Input Line Entry System)</em> adalah standar internasional untuk menuliskan struktur molekul (atom, ikatan, dan percabangannya) menjadi baris teks ASCII yang ringkas.
          </p>

          <div className="flex-grow mb-4">
            <p className="text-xs font-bold text-main uppercase tracking-wider mb-3">Simbol Dasar</p>
            <div className="space-y-2.5">
              <div className="flex items-start gap-3">
                <span className="font-mono text-xs font-bold text-blue-600 dark:text-blue-400 bg-blue-50 dark:bg-blue-900/30 rounded px-2 py-0.5 w-14 text-center flex-shrink-0">C, O, N</span>
                <p className="text-xs text-muted leading-relaxed">Huruf mewakili atom. Huruf besar = atom di luar cincin aromatik, huruf kecil (mis. <code className="font-mono">c</code>) = atom pada cincin aromatik.</p>
              </div>
              <div className="flex items-start gap-3">
                <span className="font-mono text-xs font-bold text-blue-600 dark:text-blue-400 bg-blue-50 dark:bg-blue-900/30 rounded px-2 py-0.5 w-14 text-center flex-shrink-0">( )</span>
                <p className="text-xs text-muted leading-relaxed">Tanda kurung menandai percabangan rantai dari atom sebelumnya, mis. <code className="font-mono">CC(C)C</code> untuk isobutana.</p>
              </div>
              <div className="flex items-start gap-3">
                <span className="font-mono text-xs font-bold text-blue-600 dark:text-blue-400 bg-blue-50 dark:bg-blue-900/30 rounded px-2 py-0.5 w-14 text-center flex-shrink-0">= , #</span>
                <p className="text-xs text-muted leading-relaxed">Tanda ikatan: <code className="font-mono">=</code> untuk ikatan rangkap dua, <code className="font-mono">#</code> untuk ikatan rangkap tiga. Ikatan tunggal tidak perlu simbol.</p>
              </div>
              <div className="flex items-start gap-3">
                <span className="font-mono text-xs font-bold text-blue-600 dark:text-blue-400 bg-blue-50 dark:bg-blue-900/30 rounded px-2 py-0.5 w-14 text-center flex-shrink-0">1, 2...</span>
                <p className="text-xs text-muted leading-relaxed">Angka menandai penutupan cincin — dua atom dengan angka yang sama saling terhubung, mis. <code className="font-mono">C1CCCCC1</code> untuk sikloheksana.</p>
              </div>
            </div>
          </div>

          <div className="bg-slate-100 dark:bg-slate-800 rounded-xl p-4 text-slate-900 dark:text-slate-200 font-mono text-xs border border-slate-100 dark:border-slate-700">
            <div className="flex justify-between mb-2 border-b border-slate-300 dark:border-slate-600 pb-2">
              <span>Molekul</span>
              <span>SMILES Code</span>
            </div>
            <div className="space-y-1">
              <div className="flex justify-between">
                <span className="text-black dark:text-white">Ethanol</span>
                <span className="text-green-600 dark:text-green-400 font-bold">CCO</span>
              </div>
              <div className="flex justify-between">
                <span className="text-black dark:text-white">Benzene</span>
                <span className="text-green-600 dark:text-green-400 font-bold">c1ccccc1</span>
              </div>
              <div className="flex justify-between">
                <span className="text-black dark:text-white">Carbon Dioxide</span>
                <span className="text-green-600 dark:text-green-400 font-bold">O=C=O</span>
              </div>
            </div>
          </div>
        </div>

        {/* --- CARD 3: MEMBACA HASIL PREDIKSI (Half Width) --- */}
        <div className="bg-card rounded-2xl border border-border shadow-sm p-8 hover:shadow-md transition-shadow flex flex-col">
          <div className="flex items-center gap-4 mb-4">
            <div className="w-10 h-10 rounded-lg bg-green-50 dark:bg-green-900/30 flex items-center justify-center text-green-600 dark:text-green-400">
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z" /></svg>
            </div>
            <h2 className="text-lg font-bold text-main">Membaca Hasil Prediksi</h2>
          </div>

          <p className="text-muted text-sm mb-6 leading-relaxed">
            Setiap prediksi disertai dua elemen yang membantu Anda mengevaluasi hasilnya:
          </p>

          <div className="space-y-4">
            {/* Justification */}
            <div className="flex items-start gap-3 p-3 rounded-lg border border-blue-100 dark:border-blue-900 bg-blue-50/50 dark:bg-blue-900/20">
              <div className="w-2 h-full min-h-[2.5rem] rounded-full bg-blue-500"></div>
              <div>
                <p className="text-xs font-bold text-blue-700 dark:text-blue-400 uppercase tracking-wider">Justification</p>
                <p className="text-xs text-gray-600 dark:text-gray-400 font-medium">
                  Penjelasan naratif dari AI mengenai bagaimana gugus fungsi dan struktur rantai memengaruhi sifat yang diprediksi (mis. Tg). Gunakan bagian ini untuk menilai apakah penalaran AI masuk akal secara kimiawi.
                </p>
              </div>
            </div>

            {/* Similar Compounds % Match — ranked #1, #2, #3 */}
            <div>
              <p className="text-xs font-bold text-main uppercase tracking-wider mb-2">Similar Compounds (% Match)</p>
              <p className="text-xs text-muted mb-3 leading-relaxed">
                Senyawa paling mirip secara struktural dari database ditampilkan berurutan berdasarkan ranking, bukan berdasarkan ambang batas akurasi tetap.
              </p>

              <div className="space-y-2">
                {/* Rank #1 */}
                <div className="flex items-center gap-3 p-2.5 rounded-lg border border-green-100 dark:border-green-900 bg-green-50/50 dark:bg-green-900/20">
                  <div className="w-2 h-8 rounded-full bg-green-500"></div>
                  <div className="flex-1 flex items-center justify-between">
                    <p className="text-xs font-bold text-green-700 dark:text-green-400">#1 &middot; Match Tertinggi</p>
                  </div>
                </div>

                {/* Rank #2 */}
                <div className="flex items-center gap-3 p-2.5 rounded-lg border border-blue-100 dark:border-blue-900 bg-blue-50/50 dark:bg-blue-900/20">
                  <div className="w-2 h-8 rounded-full bg-blue-500"></div>
                  <div className="flex-1 flex items-center justify-between">
                    <p className="text-xs font-bold text-blue-700 dark:text-blue-400">#2 &middot; Match Kedua</p>
                  </div>
                </div>

                {/* Rank #3 */}
                <div className="flex items-center gap-3 p-2.5 rounded-lg border border-slate-200 dark:border-slate-700 bg-slate-50/50 dark:bg-slate-800/40">
                  <div className="w-2 h-8 rounded-full bg-slate-400"></div>
                  <div className="flex-1 flex items-center justify-between">
                    <p className="text-xs font-bold text-slate-600 dark:text-slate-300">#3 &middot; Match Ketiga</p>
                  </div>
                </div>
              </div>

              <p className="text-xs text-gray-500 dark:text-gray-400 italic mt-3">
                Semakin rendah nomor rank, semakin mirip strukturnya dengan senyawa referensi. Jika ketiga match bernilai rendah, senyawa Anda kemungkinan cukup unik/novel — prediksi lebih bergantung pada reasoning AI di bagian Justification.
              </p>
            </div>
          </div>
        </div>

      </div>
    </div>
  );
}

export default WikiPage;