import { Link, useLocation, useNavigate } from 'react-router-dom';
import { type ReactNode, useState, useEffect } from 'react';
import { signOut, onAuthStateChanged, type User } from 'firebase/auth'; 
import { auth } from '../../lib/firebase';
import logoImage from '../../assets/logo.png';

interface LayoutProps {
  children: ReactNode;
}

// --- KOMPONEN HELPER ---
const MenuItem = ({ to, label, icon, active, isOpen }: { to: string, label: string, icon: ReactNode, active?: boolean, isOpen: boolean }) => (
  <Link 
    to={to} 
    title={!isOpen ? label : ''}
    className={`flex items-center px-6 py-3 text-sm font-medium transition-all duration-300 ${
      active ? 'text-blue-600 bg-blue-50 border-r-4 border-blue-600' : 'text-gray-500 hover:text-gray-900 hover:bg-gray-50'
    } ${!isOpen ? 'justify-center px-2' : ''}`}
  >
    <div className="w-6 h-6 flex-shrink-0">{icon}</div>
    <span className={`ml-4 transition-opacity duration-300 whitespace-nowrap ${isOpen ? 'opacity-100' : 'opacity-0 w-0 hidden'}`}>{label}</span>
  </Link>
);

const MenuButton = ({ label, icon, onClick, isOpen }: { label: string, icon: ReactNode, onClick?: () => void, isOpen: boolean }) => (
  <button 
    onClick={onClick}
    title={!isOpen ? label : ''}
    className={`w-full flex items-center px-6 py-3 text-sm font-medium text-gray-500 hover:text-gray-900 hover:bg-gray-50 transition-all duration-300 text-left ${!isOpen ? 'justify-center px-2' : ''}`}
  >
    <div className="w-6 h-6 flex-shrink-0">{icon}</div>
    <span className={`ml-4 transition-opacity duration-300 whitespace-nowrap ${isOpen ? 'opacity-100' : 'opacity-0 w-0 hidden'}`}>{label}</span>
  </button>
);

// --- KOMPONEN UTAMA ---
function DashboardLayout({ children }: LayoutProps) {
  const location = useLocation();
  const navigate = useNavigate();
  const isActive = (path: string) => location.pathname === path;
  const [isSidebarOpen, setIsSidebarOpen] = useState(true);
  
  // State untuk User
  const [user, setUser] = useState<User | null>(null);

  // 1. CEK STATUS LOGIN
  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (currentUser) => {
      setUser(currentUser);
    });
    return () => unsubscribe();
  }, []);

  // 2. FUNGSI LOGOUT
  const handleLogout = async () => {
    try {
      await signOut(auth);
      navigate('/login');
    } catch (error) {
      console.error("Gagal logout:", error);
    }
  };

  return (
    <div className="flex h-screen bg-gray-50 overflow-hidden">
      <aside className={`${isSidebarOpen ? 'w-64' : 'w-20'} bg-white border-r border-gray-200 hidden md:flex flex-col justify-between transition-all duration-300 ease-in-out`}>
        
        {/* Header Logo */}
        <div className={`h-16 flex items-center border-b border-gray-200 shrink-0 transition-all duration-300 ${isSidebarOpen ? 'px-6 justify-between' : 'px-0 justify-center'}`}>
            <div className="flex items-center overflow-hidden">
              <div className="h-10 w-10 rounded-full overflow-hidden border border-gray-200 flex items-center justify-center bg-white flex-shrink-0">
                  <img src={logoImage} alt="Logo" className="h-full w-full object-cover" />
              </div>
              <span className={`font-bold text-gray-900 text-lg tracking-tight ml-3 transition-opacity duration-300 whitespace-nowrap ${isSidebarOpen ? 'opacity-100' : 'opacity-0 w-0 hidden'}`}>
                NovaChem
              </span>
            </div>
             <button onClick={() => setIsSidebarOpen(!isSidebarOpen)} className={`text-gray-400 hover:text-gray-600 transition-colors ${!isSidebarOpen ? 'hidden' : 'block'}`}>
               <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11 19l-7-7 7-7m8 14l-7-7 7-7" /></svg>
             </button>
        </div>
        
        {!isSidebarOpen && (
           <button onClick={() => setIsSidebarOpen(true)} className="w-full flex justify-center py-2 text-gray-400 hover:text-blue-600 border-b border-gray-100">
              <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 5l7 7-7 7M5 5l7 7-7 7" /></svg>
           </button>
        )}

        {/* Menu Items */}
        <div className="flex-1 overflow-y-auto py-4 overflow-x-hidden">
          
          <div className="space-y-1">
             {/* HOME: Rumah */}
             <MenuItem isOpen={isSidebarOpen} to="/" label="Home" active={isActive('/')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M3 12l2-2m0 0l7-7 7 7M5 10v10a1 1 0 001 1h3m10-11l2 2m-2-2v10a1 1 0 01-1 1h-3m-6 0a1 1 0 001-1v-4a1 1 0 011-1h2a1 1 0 011 1v4a1 1 0 001 1m-6 0h6" /></svg>} />
             
             {/* MODULES: Grid 4 Kotak (Beda) */}
             <MenuItem isOpen={isSidebarOpen} to="/explore" label="Modules" active={isActive('/explore')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M4 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2V6zM14 6a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2V6zM4 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2H6a2 2 0 01-2-2v-2zM14 16a2 2 0 012-2h2a2 2 0 012 2v2a2 2 0 01-2 2h-2a2 2 0 01-2-2v-2z" /></svg>} />
             
             {/* LAB TOOLS: Tabung Reaksi (Spesifik) */}
             <MenuItem isOpen={isSidebarOpen} to="#" label="Lab Tools" icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" /></svg>} />
          </div>

          <div className="my-4 border-t border-gray-200 mx-4"></div>

          <div className="space-y-1">
            {/* LIBRARY: Bookmark */}
            <MenuItem isOpen={isSidebarOpen} to="/library" label="Library" active={isActive('/library')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" /></svg>} />
            
            {/* WIKI: Buku Terbuka */}
            <MenuItem isOpen={isSidebarOpen} to="/wiki" label="Wiki's" active={isActive('/wiki')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" /></svg>} />
            
            {/* ABOUT: Info Circle (Tetap) */}
            <MenuItem isOpen={isSidebarOpen} to="/about" label="About" active={isActive('/about')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>} />

            {/* HISTORY: Jam (Time) */}
            <MenuItem isOpen={isSidebarOpen} to="/history" label="History" active={isActive('/history')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>} />
            
            {/* SEARCH: Kaca Pembesar (Tetap) */}
            <MenuItem isOpen={isSidebarOpen} to="#" label="Search" icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" /></svg>} />
          </div>
          
          <div className="my-4 border-t border-gray-200 mx-4"></div>

          <div className="space-y-1">
            <MenuButton isOpen={isSidebarOpen} label="Dark Mode" icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M20.354 15.354A9 9 0 018.646 3.646 9.003 9.003 0 0012 21a9.003 9.003 0 008.354-5.646z" /></svg>} />
          </div>
        </div>

        {/* FOOTER USER & LOGOUT */}
        <div className="p-4 border-t border-gray-200 space-y-1">
           <MenuItem 
             isOpen={isSidebarOpen} 
             to="#" 
             label={user ? (user.displayName || user.email || 'User') : 'Guest'} 
             icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z" /></svg>} 
           />
           <MenuButton 
             isOpen={isSidebarOpen} 
             label={user ? "Logout" : "Login"} 
             onClick={user ? handleLogout : () => navigate('/login')} 
             icon={user ? <svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M17 16l4-4m0 0l-4-4m4 4H7m6 4v1a3 3 0 01-3 3H6a3 3 0 01-3-3V7a3 3 0 013-3h4a3 3 0 013 3v1" /></svg> : <svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M11 16l-4-4m0 0l4-4m-4 4h14m-5 4v1a3 3 0 01-3 3H6a3 3 0 01-3-3V7a3 3 0 013-3h7a3 3 0 013 3v1" /></svg>} 
           />
        </div>
      </aside>

      {/* CONTENT */}
      <main className="flex-1 overflow-y-auto relative">
        <div className="md:hidden absolute top-4 left-4 z-20">
           <button className="p-2 bg-white rounded-md shadow-md border border-gray-200 text-gray-600">☰</button>
        </div>
        <div className="p-8 min-h-full">{children}</div>
      </main>
    </div>
  );
}

export default DashboardLayout;