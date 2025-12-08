function WikiPage() {
  return (
    <div className="max-w-4xl mx-auto">
      <h1 className="text-3xl font-bold text-gray-900 mb-8">Dokumentasi & Panduan</h1>

      <div className="grid grid-cols-1 gap-8">
        
        {/* Topik 1: Cara Menggunakan */}
        <div className="bg-white p-8 rounded-xl shadow-sm border border-gray-200">
          <h2 className="text-xl font-bold text-blue-600 mb-4">🚀 Cara Memulai Pencarian</h2>
          <p className="text-gray-600 mb-4">
            NovaChemAI memungkinkan Anda mencari prediksi properti polimer dengan dua cara:
          </p>
          <ol className="list-decimal list-inside space-y-2 text-gray-700 ml-2">
            <li><strong>Nama Umum:</strong> Masukkan nama polimer seperti "Polyethylene" atau "PET".</li>
            <li><strong>Format SMILES:</strong> Masukkan notasi kimia linear untuk struktur spesifik.</li>
          </ol>
        </div>

        {/* Topik 2: Apa itu SMILES */}
        <div className="bg-white p-8 rounded-xl shadow-sm border border-gray-200">
          <h2 className="text-xl font-bold text-blue-600 mb-4">🧪 Apa itu Format SMILES?</h2>
          <p className="text-gray-600 mb-4">
            <strong>SMILES</strong> (Simplified Molecular Input Line Entry System) adalah cara merepresentasikan struktur kimia menggunakan karakter teks ASCII.
          </p>
          <div className="bg-slate-100 p-4 rounded-lg border border-slate-200 font-mono text-sm mb-4">
            Contoh: CCO (Ethanol), c1ccccc1 (Benzene)
          </div>
          <p className="text-gray-600">
            Sistem AI kami membaca format ini untuk memahami ikatan atom dan memprediksi sifat material.
          </p>
        </div>

        {/* Topik 3: Interpretasi Skor */}
        <div className="bg-white p-8 rounded-xl shadow-sm border border-gray-200">
          <h2 className="text-xl font-bold text-blue-600 mb-4">📊 Menginterpretasikan Skor AI</h2>
          <p className="text-gray-600 mb-4">
            Setiap prediksi dilengkapi dengan <strong>Confidence Score</strong> (0% - 100%).
          </p>
          <ul className="space-y-2 text-gray-700">
            <li><span className="font-bold text-green-600">90% - 100%:</span> Sangat Akurat (Data validasi tersedia).</li>
            <li><span className="font-bold text-yellow-600">70% - 89%:</span> Prediksi Moderat (Perlu verifikasi lab).</li>
            <li><span className="font-bold text-red-600">&lt; 70%:</span> Prediksi Rendah (Struktur tidak umum).</li>
          </ul>
        </div>

      </div>
    </div>
  );
}

export default WikiPage;