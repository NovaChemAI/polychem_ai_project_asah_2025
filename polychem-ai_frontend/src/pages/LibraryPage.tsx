import { useEffect, useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { auth } from '../lib/firebase';
import { onAuthStateChanged } from 'firebase/auth';
import { getUserLibrary, removeFromLibrary, type SavedChemical } from '../services/dbService';

// Define your backend URL constant
const BACKEND_URL = "https://slim-danika-polychem-ab276767.koyeb.app";

function LibraryPage() {
  const [items, setItems] = useState<SavedChemical[]>([]);
  const [loading, setLoading] = useState(true);
  const [userUid, setUserUid] = useState<string | null>(null);
  const navigate = useNavigate();

  // --- STATE FOR DELETE MODAL ---
  const [isModalOpen, setIsModalOpen] = useState(false);
  const [itemToDelete, setItemToDelete] = useState<string | null>(null);
  const [isDeleting, setIsDeleting] = useState(false);

  // --- 1. Fetch Data Function ---
  const fetchItems = async (uid: string) => {
    // Note: Jika ingin skeleton muncul setiap kali fetch, uncomment baris bawah
    // setLoading(true); 
    const data = await getUserLibrary(uid);
    setItems(data);
    setLoading(false);
  };

  // --- 2. Check Login & Load Data ---
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, async (user) => {
      if (user) {
        setUserUid(user.uid);
        fetchItems(user.uid);
      } else {
        setItems([]);
        setLoading(false);
      }
    });
    return () => unsubscribe();
  }, [navigate]);

  // --- 3. Trigger Modal (Open) ---
  const handleDeleteClick = (itemId: string, e: React.MouseEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setItemToDelete(itemId);
    setIsModalOpen(true);
  };

  // --- 4. Confirm Delete Action ---
  const confirmDelete = async () => {
    if (!userUid || !itemToDelete) return;

    setIsDeleting(true);
    try {
      setItems(prev => prev.filter(item => item.id !== itemToDelete));
      const success = await removeFromLibrary(userUid, itemToDelete);
      
      if (!success) {
        alert("Gagal menghapus data dari server. Data akan dimuat ulang.");
        fetchItems(userUid);
      }
    } catch (error) {
      console.error(error);
      alert("Terjadi kesalahan sistem.");
    } finally {
      setIsDeleting(false);
      setIsModalOpen(false);
      setItemToDelete(null);
    }
  };

  // --- 5. Cancel Delete ---
  const cancelDelete = () => {
    setIsModalOpen(false);
    setItemToDelete(null);
  };

  // Helper function to get full image URL
  const getImageUrl = (imagePath: string | undefined) => {
    if (!imagePath) return null;
    if (imagePath.startsWith('http')) return imagePath; 
    return `${BACKEND_URL}${imagePath.startsWith('/') ? '' : '/'}${imagePath}`;
  };


  // --- UI: LOADING SKELETON (DIPERBAIKI) ---
  // Menggunakan Grid Skeleton, bukan Spinner bulat
  if (loading) {
    return (
      <div className="max-w-7xl mx-auto px-4 py-8 animate-pulse">
        {/* Header Skeleton */}
        <div className="mb-8 flex justify-between items-end">
           <div>
              <div className="h-8 bg-gray-200 dark:bg-slate-700 w-48 rounded mb-2"></div>
              <div className="h-4 bg-gray-200 dark:bg-slate-700 w-64 rounded"></div>
           </div>
           <div className="h-8 bg-gray-200 dark:bg-slate-700 w-20 rounded-lg"></div>
        </div>

        {/* Cards Grid Skeleton */}
        <div className="grid grid-cols-1 gap-6 sm:grid-cols-2 lg:grid-cols-3">
           {/* Render 6 kartu palsu */}
           {[1, 2, 3, 4, 5, 6].map((i) => (
             <div key={i} className="bg-white dark:bg-slate-800 rounded-xl border border-gray-200 dark:border-slate-700 p-6 h-64 flex flex-col justify-between">
                {/* Top: Image Placeholder & Button */}
                <div className="flex justify-between items-start mb-4">
                   <div className="w-16 h-16 bg-gray-200 dark:bg-slate-700 rounded-lg"></div>
                   <div className="w-8 h-8 bg-gray-200 dark:bg-slate-700 rounded-full"></div>
                </div>

                {/* Middle: Text Lines */}
                <div>
                  <div className="h-6 bg-gray-200 dark:bg-slate-700 rounded w-3/4 mb-3"></div>
                  <div className="h-8 bg-gray-200 dark:bg-slate-700 rounded w-full mb-2"></div>
                </div>

                {/* Bottom: Footer Info */}
                <div className="flex justify-between pt-3 border-t border-gray-100 dark:border-slate-700">
                   <div className="h-4 bg-gray-200 dark:bg-slate-700 w-1/3 rounded"></div>
                   <div className="h-4 bg-gray-200 dark:bg-slate-700 w-1/4 rounded"></div>
                </div>
             </div>
           ))}
        </div>
      </div>
    );
  }

  // --- Main UI ---
  return (
    <div className="max-w-7xl mx-auto px-4 py-8 relative">
      
      {/* --- MODAL CONFIRMATION --- */}
      {isModalOpen && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/50 backdrop-blur-sm p-4">
          <div className="bg-white dark:bg-slate-800 rounded-2xl shadow-2xl w-full max-w-sm border border-gray-200 dark:border-slate-700 transform transition-all scale-100 p-6">
            <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
              Hapus Item?
            </h3>
            <p className="text-gray-500 dark:text-slate-400 mb-6">
              Tindakan ini tidak dapat dibatalkan. Item akan dihapus permanen dari library Anda.
            </p>
            <div className="flex justify-end gap-3">
              <button
                onClick={cancelDelete}
                disabled={isDeleting}
                className="px-4 py-2 rounded-lg text-gray-700 dark:text-slate-300 hover:bg-gray-100 dark:hover:bg-slate-700 font-medium transition-colors"
              >
                Batal
              </button>
              <button
                onClick={confirmDelete}
                disabled={isDeleting}
                className="px-4 py-2 rounded-lg bg-red-600 hover:bg-red-700 text-white font-bold shadow-lg shadow-red-500/30 transition-all active:scale-95 disabled:opacity-70 flex items-center gap-2"
              >
                {isDeleting ? (
                  <>
                    <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin"></div>
                    Menghapus...
                  </>
                ) : (
                  "Ya, Hapus"
                )}
              </button>
            </div>
          </div>
        </div>
      )}

      <div className="mb-8 flex justify-between items-end">
        <div>
           <h1 className="text-2xl font-bold text-gray-900 dark:text-white">My Library</h1>
           <p className="text-gray-600 dark:text-slate-400">Koleksi senyawa dan hasil AI yang tersimpan</p>
        </div>
        <div className="bg-blue-50 dark:bg-blue-900/30 text-blue-700 dark:text-blue-400 px-3 py-1 rounded-lg text-sm font-bold border border-blue-100 dark:border-blue-800">
            {items.length} Items
        </div>
      </div>

      {items.length > 0 ? (
        <div className="grid grid-cols-1 gap-6 sm:grid-cols-2 lg:grid-cols-3">
          {items.map((item) => {
             const imageSrc = getImageUrl(item.image_url || item.image);

             return (
            <div 
              key={item.id} 
              className="bg-white dark:bg-slate-800 block overflow-hidden shadow-sm rounded-xl border border-gray-200 dark:border-slate-700 hover:border-blue-400 dark:hover:border-blue-500 hover:shadow-md transition-all relative group"
            >
              <div className="p-6">
                <div className="flex justify-between items-start mb-4">
                  <div className="relative">
                      {imageSrc ? (
                        <img 
                            src={imageSrc} 
                            alt={item.name} 
                            className="w-16 h-16 object-contain border dark:border-slate-600 bg-white rounded p-1" 
                            onError={(e) => {
                                (e.target as HTMLImageElement).style.display = 'none';
                                (e.target as HTMLImageElement).nextElementSibling?.classList.remove('hidden');
                            }}
                        />
                      ) : (
                        <div className="inline-flex items-center justify-center w-16 h-16 rounded-lg bg-blue-50 dark:bg-blue-900/50 text-blue-700 dark:text-blue-300 font-bold text-2xl border border-blue-100 dark:border-blue-800">
                            {item.name ? item.name.charAt(0).toUpperCase() : '?'}
                        </div>
                      )}
                      
                      <div className="hidden inline-flex items-center justify-center w-16 h-16 rounded-lg bg-blue-50 dark:bg-blue-900/50 text-blue-700 dark:text-blue-300 font-bold text-2xl absolute top-0 left-0 border border-blue-100 dark:border-blue-800">
                           {item.name ? item.name.charAt(0).toUpperCase() : '?'}
                      </div>
                      
                      {item.isAiResult && (
                        <span className="absolute -top-2 -right-2 bg-purple-100 dark:bg-purple-900/60 text-purple-700 dark:text-purple-300 text-[10px] font-bold px-1.5 py-0.5 rounded-full border border-purple-200 dark:border-purple-700 shadow-sm">
                            AI
                        </span>
                      )}
                  </div>
                  
                  {/* DELETE BUTTON */}
                  <button 
                    onClick={(e) => handleDeleteClick(item.id, e)}
                    className="p-2 text-gray-400 dark:text-slate-500 hover:text-red-600 dark:hover:text-red-400 hover:bg-red-50 dark:hover:bg-red-900/20 rounded-full transition-colors z-10"
                    title="Hapus"
                  >
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" /></svg>
                  </button>
                </div>
                
                <h3 className="text-lg font-bold text-gray-900 dark:text-white mb-1 truncate" title={item.name}>
                    {item.name}
                </h3>
                
                <p className="text-xs text-gray-500 dark:text-slate-400 font-mono bg-gray-50 dark:bg-slate-900/50 border border-gray-100 dark:border-slate-700 p-2 rounded truncate mb-3" title={item.smiles}>
                  {item.smiles}
                </p>

                <div className="flex justify-between items-center text-xs text-gray-400 dark:text-slate-500 border-t border-gray-100 dark:border-slate-700 pt-3">
                  <span>
                    {item.savedAt?.seconds 
                        ? new Date(item.savedAt.seconds * 1000).toLocaleDateString() 
                        : "Baru saja"}
                  </span>
                  <span className="font-semibold text-blue-600 dark:text-blue-400 bg-blue-50 dark:bg-blue-900/20 border border-blue-100 dark:border-blue-800 px-2 py-0.5 rounded">
                    {item.category || "Unknown"}
                  </span>
                </div>
              </div>
            </div>
          )})}
        </div>
      ) : (
        <div className="p-12 text-center bg-white dark:bg-slate-800 rounded-xl border border-dashed border-gray-300 dark:border-slate-700">
          <p className="text-gray-500 dark:text-slate-400 mb-4">Belum ada data di library.</p>
          <Link to="/" className="text-blue-600 dark:text-blue-400 hover:underline font-semibold">
            Mulai Mencari Senyawa
          </Link>
        </div>
      )}
    </div>
  );
}

export default LibraryPage;