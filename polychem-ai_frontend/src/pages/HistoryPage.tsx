import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { auth } from '../lib/firebase';
import { onAuthStateChanged } from 'firebase/auth';
import { getUserHistory, type HistoryItem } from '../services/dbService';

function HistoryPage() {
  const [historyItems, setHistoryItems] = useState<HistoryItem[]>([]);
  const [loading, setLoading] = useState(true);
  const navigate = useNavigate();

  useEffect(() => {
    // Simulasi delay sedikit agar skeleton terlihat (opsional, bisa dihapus timeout-nya)
    const unsubscribe = onAuthStateChanged(auth, async (user) => {
      if (user) {
        const data = await getUserHistory(user.uid);
        setHistoryItems(data);
      } else {
        navigate('/login');
      }
      setLoading(false);
    });
    return () => unsubscribe();
  }, [navigate]);

  const handleItemClick = (item: HistoryItem) => {
    try {
        const predictionData = JSON.parse(item.full_data);
        navigate('/result', { state: { predictionData: predictionData } });
    } catch (e) {
        console.error("Gagal membuka history:", e);
        alert("Data history rusak atau tidak valid.");
    }
  };

  // --- TAMPILAN SKELETON LOADING (BARU) ---
  if (loading) {
    return (
      <div className="max-w-4xl mx-auto px-4 py-8 animate-pulse">
        {/* Header Skeleton */}
        <div className="mb-6">
           <div className="h-8 bg-gray-200 dark:bg-slate-700 w-48 rounded-lg mb-2"></div>
           <div className="h-4 bg-gray-200 dark:bg-slate-700 w-64 rounded-lg"></div>
        </div>

        {/* List Skeleton */}
        <div className="bg-card rounded-xl border border-border overflow-hidden">
           <div className="divide-y divide-border">
              {[1, 2, 3, 4, 5].map((i) => (
                <div key={i} className="p-4 flex justify-between items-center">
                   <div className="flex items-center gap-4">
                      <div className="w-10 h-10 bg-gray-200 dark:bg-slate-700 rounded-full"></div>
                      <div className="space-y-2">
                         <div className="h-4 bg-gray-200 dark:bg-slate-700 w-32 rounded"></div>
                         <div className="h-3 bg-gray-200 dark:bg-slate-700 w-48 rounded"></div>
                      </div>
                   </div>
                   <div className="w-6 h-6 bg-gray-200 dark:bg-slate-700 rounded"></div>
                </div>
              ))}
           </div>
        </div>
      </div>
    );
  }

  // --- TAMPILAN UTAMA ---
  return (
    <div className="max-w-4xl mx-auto px-4 py-8">
      <div className="mb-6">
        {/* Update: text-main */}
        <h1 className="text-2xl font-bold text-main">Research History</h1>
        {/* Update: text-muted */}
        <p className="text-muted">Jejak analisis otomatis Anda.</p>
      </div>

      {/* Update: bg-card, border-border */}
      <div className="bg-card rounded-xl shadow-sm border border-border overflow-hidden">
        {historyItems.length > 0 ? (
          // Update: divide-border
          <ul className="divide-y divide-border">
            {historyItems.map((item) => (
              <li 
                key={item.id} 
                onClick={() => handleItemClick(item)}
                // Update: hover states for dark mode
                className="hover:bg-blue-50 dark:hover:bg-blue-900/10 transition-colors cursor-pointer p-4 flex justify-between items-center group"
              >
                <div className="flex items-center gap-4">
                    {/* Update: Icon background colors */}
                    <div className="w-10 h-10 bg-gray-100 dark:bg-slate-800 rounded-full flex items-center justify-center text-gray-500 dark:text-slate-400 group-hover:bg-blue-200 dark:group-hover:bg-blue-900/40 group-hover:text-blue-700 dark:group-hover:text-blue-400 transition-colors">
                        <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={1.5} stroke="currentColor" className="w-5 h-5">
                            <path strokeLinecap="round" strokeLinejoin="round" d="M12 6v6h4.5m4.5 0a9 9 0 11-18 0 9 9 0 0118 0z" />
                        </svg>
                    </div>
                    <div>
                        {/* Update: text-main */}
                        <p className="font-bold text-main">{item.result_name}</p>
                        {/* Update: text-muted */}
                        <p className="text-xs text-muted font-mono">Input: {item.query}</p>
                        <p className="text-[10px] text-gray-400 dark:text-slate-500 mt-1">
                           {item.timestamp?.seconds 
                             ? new Date(item.timestamp.seconds * 1000).toLocaleString('id-ID', { dateStyle: 'medium', timeStyle: 'short' }) 
                             : "Baru saja"}
                        </p>
                    </div>
                </div>
                <div className="flex items-center gap-2">
                    {/* Update: text-blue-400 for dark mode */}
                    <span className="text-sm text-blue-600 dark:text-blue-400 font-medium opacity-0 group-hover:opacity-100 transition-opacity">
                        Lihat Hasil
                    </span>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={2} stroke="currentColor" className="w-4 h-4 text-gray-300 dark:text-slate-600 group-hover:text-blue-600 dark:group-hover:text-blue-400">
                        <path strokeLinecap="round" strokeLinejoin="round" d="M8.25 4.5l7.5 7.5-7.5 7.5" />
                    </svg>
                </div>
              </li>
            ))}
          </ul>
        ) : (
          <div className="p-12 text-center text-muted">
             <p className="mb-4">Belum ada riwayat pencarian.</p>
          </div>
        )}
      </div>
    </div>
  );
}

export default HistoryPage;