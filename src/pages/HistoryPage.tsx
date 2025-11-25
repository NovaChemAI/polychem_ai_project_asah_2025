import { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import { auth } from '../lib/firebase';
import { getUserHistory } from '../services/dbService';
import { onAuthStateChanged } from 'firebase/auth';

interface HistoryItem {
  keyword: string;
  timestamp: string;
}

function HistoryPage() {
  const [history, setHistory] = useState<HistoryItem[]>([]);
  const [loading, setLoading] = useState(true);
  const [userName, setUserName] = useState('');

  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, async (user) => {
      if (user) {
        setUserName(user.displayName || 'User');
        const data = await getUserHistory(user.uid);
        setHistory(data);
      }
      setLoading(false);
    });

    return () => unsubscribe();
  }, []);

  return (
    <div className="max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-2xl font-bold text-gray-900">Riwayat Pencarian</h1>
        <p className="text-gray-600">Daftar aktivitas pencarian senyawa oleh {userName}.</p>
      </div>

      <div className="bg-white rounded-xl shadow-sm border border-gray-200 overflow-hidden">
        {loading ? (
          <div className="p-8 text-center text-gray-500">Memuat data...</div>
        ) : history.length > 0 ? (
          <ul className="divide-y divide-gray-100">
            {history.map((item, index) => (
              <li key={index} className="hover:bg-gray-50 transition-colors">
                <Link 
                  // Arahkan kembali ke halaman explore saat diklik
                  to="/explore" 
                  className="flex items-center justify-between px-6 py-4"
                >
                  <div className="flex items-center gap-4">
                    <div className="w-10 h-10 rounded-full bg-blue-50 flex items-center justify-center text-blue-600">
                      {/* Ikon Jam */}
                      <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>
                    </div>
                    <div>
                      <p className="font-medium text-gray-900">{item.keyword}</p>
                      <p className="text-xs text-gray-400">
                        {/* Format Tanggal Sederhana */}
                        {new Date(item.timestamp).toLocaleDateString('id-ID', { 
                          day: 'numeric', month: 'long', hour: '2-digit', minute: '2-digit' 
                        })}
                      </p>
                    </div>
                  </div>
                  <span className="text-blue-600 text-sm font-medium">Lihat Hasil &rarr;</span>
                </Link>
              </li>
            ))}
          </ul>
        ) : (
          // Tampilan jika History Kosong
          <div className="p-12 text-center flex flex-col items-center">
            <div className="w-16 h-16 bg-gray-100 rounded-full flex items-center justify-center mb-4 text-gray-400">
              <svg className="w-8 h-8" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" /></svg>
            </div>
            <h3 className="text-lg font-medium text-gray-900">Belum ada riwayat</h3>
            <p className="text-gray-500 mt-1 mb-6">Mulai eksplorasi senyawa untuk melihat riwayat Anda di sini.</p>
            <Link to="/explore" className="bg-gray-900 text-white px-6 py-2 rounded-lg text-sm font-medium hover:bg-gray-800 transition-colors">
              Mulai Explorasi
            </Link>
          </div>
        )}
      </div>
    </div>
  );
}

export default HistoryPage;