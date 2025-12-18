import { Link, useLocation, useNavigate } from 'react-router-dom';
import { type ReactNode, useState, useEffect } from 'react';
import { signOut, onAuthStateChanged, type User } from 'firebase/auth';
import { auth } from '../../lib/firebase';
import logoImage from '../../assets/logo.png';
import { useTheme } from '../../context/ThemeContext';
import ConfirmationModal from '../ui/ConfirmationModal';
import toast from 'react-hot-toast';

interface LayoutProps {
  children: ReactNode;
}

// Komponen Navigasi Reusable dengan ICON ASLI
const NavLinks = ({ isOpen, isActive, onItemClick }: { isOpen: boolean, isActive: (path: string) => boolean, onItemClick?: () => void }) => (
  <div className="space-y-1">
    <MenuItem onItemClick={onItemClick} isOpen={isOpen} to="/" label="Home" active={isActive('/')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M3 12l2-2m0 0l7-7 7 7M5 10v10a1 1 0 001 1h3m10-11l2 2m-2-2v10a1 1 0 01-1 1h-3m-6 0a1 1 0 001-1v-4a1 1 0 011-1h2a1 1 0 011 1v4a1 1 0 001 1m-6 0h6" /></svg>} />
    <MenuItem onItemClick={onItemClick} isOpen={isOpen} to="/library" label="Library" active={isActive('/library')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M5 5a2 2 0 012-2h10a2 2 0 012 2v16l-7-3.5L5 21V5z" /></svg>} />
    <MenuItem onItemClick={onItemClick} isOpen={isOpen} to="/wiki" label="Wiki's" active={isActive('/wiki')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" /></svg>} />
    <MenuItem onItemClick={onItemClick} isOpen={isOpen} to="/about" label="About" active={isActive('/about')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>} />
    <MenuItem onItemClick={onItemClick} isOpen={isOpen} to="/history" label="History" active={isActive('/history')} icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>} />
  </div>
);

const MenuItem = ({ to, label, icon, active, isOpen, onItemClick }: { to: string, label: string, icon: ReactNode, active?: boolean, isOpen: boolean, onItemClick?: () => void }) => (
  <Link
    to={to}
    onClick={onItemClick}
    className={`flex items-center px-6 py-3 text-sm font-medium transition-all duration-300
    ${active 
      ? 'text-blue-600 bg-blue-50 border-r-4 border-blue-600 dark:bg-blue-900/20 dark:text-blue-400' 
      : 'text-gray-500 hover:text-gray-900 hover:bg-gray-50 dark:text-gray-400 dark:hover:text-gray-100 dark:hover:bg-slate-700'
    } ${!isOpen ? 'justify-center px-2' : ''}`}
  >
    <div className="w-6 h-6 flex-shrink-0">{icon}</div>
    <span className={`ml-4 transition-opacity duration-300 whitespace-nowrap ${isOpen ? 'opacity-100' : 'opacity-0 w-0 hidden'}`}>{label}</span>
  </Link>
);

const MenuButton = ({ label, icon, onClick, isOpen }: { label: string, icon: ReactNode, onClick?: () => void, isOpen: boolean }) => (
  <button
    onClick={onClick}
    className={`w-full flex items-center px-6 py-3 text-sm font-medium text-gray-500 hover:text-gray-900 hover:bg-gray-50 dark:text-gray-400 dark:hover:text-gray-100 dark:hover:bg-slate-700 transition-all duration-300 text-left ${!isOpen ? 'justify-center px-2' : ''}`}
  >
    <div className="w-6 h-6 flex-shrink-0">{icon}</div>
    <span className={`ml-4 transition-opacity duration-300 whitespace-nowrap ${isOpen ? 'opacity-100' : 'opacity-0 w-0 hidden'}`}>{label}</span>
  </button>
);

function DashboardLayout({ children }: LayoutProps) {
  const location = useLocation();
  const navigate = useNavigate();
  const isActive = (path: string) => location.pathname === path;
  
  const [isSidebarOpen, setIsSidebarOpen] = useState(true);
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);
  const { theme, toggleTheme } = useTheme();
  const [user, setUser] = useState<User | null>(null);
  const [isLogoutModalOpen, setIsLogoutModalOpen] = useState(false);
  const [isLoggingOut, setIsLoggingOut] = useState(false);

  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (currentUser) => setUser(currentUser));
    return () => unsubscribe();
  }, []);

  useEffect(() => {
    setIsMobileMenuOpen(false);
  }, [location.pathname]);

  const handleLogoutClick = () => setIsLogoutModalOpen(true);

  const handleConfirmLogout = async () => {
    setIsLoggingOut(true);
    try {
      await signOut(auth);
      toast.success("Berhasil Logout");
      navigate('/login');
    } catch { 
      // Catch tanpa variabel untuk menghindari error ESLint unused-vars
      toast.error("Gagal Logout");
    } finally {
      setIsLoggingOut(false);
      setIsLogoutModalOpen(false);
    }
  };

  return (
    <div className="flex h-screen bg-gray-50 dark:bg-slate-900 overflow-hidden transition-colors duration-300">
      
      {/* SIDEBAR DESKTOP */}
      <aside className={`${isSidebarOpen ? 'w-64' : 'w-20'} bg-white dark:bg-slate-800 border-r border-gray-200 dark:border-slate-700 hidden md:flex flex-col justify-between transition-all duration-300 ease-in-out`}>
        <div className={`h-16 flex items-center border-b border-gray-200 dark:border-slate-700 shrink-0 transition-all duration-300 ${isSidebarOpen ? 'px-6 justify-between' : 'px-0 justify-center'}`}>
            <div className="flex items-center overflow-hidden">
              <img src={logoImage} alt="Logo" className="h-10 w-10 object-contain" />
              <span className={`font-bold text-gray-900 dark:text-white text-lg tracking-tight ml-3 transition-opacity duration-300 whitespace-nowrap ${isSidebarOpen ? 'opacity-100' : 'opacity-0 w-0 hidden'}`}>PolyChemAI</span>
            </div>
            {isSidebarOpen && (
              <button onClick={() => setIsSidebarOpen(false)} className="text-gray-400 hover:text-gray-600 dark:hover:text-gray-200 transition-colors">
                <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11 19l-7-7 7-7m8 14l-7-7 7-7" /></svg>
              </button>
            )}
        </div>

        <div className="flex-1 overflow-y-auto py-4 overflow-x-hidden">
          <NavLinks isOpen={isSidebarOpen} isActive={isActive} />
          <div className="my-4 border-t border-gray-200 dark:border-slate-700 mx-4"></div>
          <MenuButton 
            isOpen={isSidebarOpen} 
            label={theme === 'dark' ? 'Light Mode' : 'Dark Mode'} 
            onClick={toggleTheme} 
            icon={theme === 'dark' ? (
              <svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24" className="text-yellow-400"><path strokeLinecap="round" strokeLinejoin="round" d="M12 3v1m0 16v1m9-9h-1M4 12H3m15.364 6.364l-.707-.707M6.343 6.343l-.707-.707m12.728 0l-.707.707M6.343 17.657l-.707.707M16 12a4 4 0 11-8 0 4 4 0 018 0z" /></svg>
            ) : (
              <svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M20.354 15.354A9 9 0 018.646 3.646 9.003 9.003 0 0012 21a9.003 9.003 0 008.354-5.646z" /></svg>
            )} 
          />
        </div>

        <div className="p-4 border-t border-gray-200 dark:border-slate-700 space-y-1">
           <MenuItem
             isOpen={isSidebarOpen}
             to="#"
             label={user ? (user.displayName || user.email || 'User') : 'Guest'}
             icon={<svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z" /></svg>}
           />
           <MenuButton 
             isOpen={isSidebarOpen} 
             label={user ? "Logout" : "Login"} 
             onClick={user ? handleLogoutClick : () => navigate('/login')} 
             icon={user ? <svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M17 16l4-4m0 0l-4-4m4 4H7m6 4v1a3 3 0 01-3 3H6a3 3 0 01-3-3V7a3 3 0 013-3h4a3 3 0 013 3v1" /></svg> : <svg fill="none" stroke="currentColor" strokeWidth="2" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" d="M11 16l-4-4m0 0l4-4m-4 4h14m-5 4v1a3 3 0 01-3 3H6a3 3 0 01-3-3V7a3 3 0 013-3h7a3 3 0 013 3v1" /></svg>} 
           />
        </div>
      </aside>

      {/* MOBILE DRAWER */}
      {isMobileMenuOpen && (
        <div className="fixed inset-0 z-40 md:hidden">
          <div className="fixed inset-0 bg-black/50 backdrop-blur-sm" onClick={() => setIsMobileMenuOpen(false)}></div>
          <nav className="fixed top-0 left-0 bottom-0 w-64 bg-white dark:bg-slate-800 shadow-xl z-50 flex flex-col">
            <div className="h-16 flex items-center px-6 border-b border-gray-100 dark:border-slate-700">
              <span className="font-bold text-blue-600 text-xl">PolyChemAI</span>
            </div>
            <div className="flex-1 py-4 overflow-y-auto">
              <NavLinks isOpen={true} isActive={isActive} onItemClick={() => setIsMobileMenuOpen(false)} />
            </div>
          </nav>
        </div>
      )}

      {/* MAIN CONTENT */}
      <main className="flex-1 overflow-y-auto relative dark:text-gray-100">
        <div className="md:hidden flex items-center justify-between p-4 bg-white dark:bg-slate-800 border-b border-gray-200 dark:border-slate-700 sticky top-0 z-20">
            <div className="flex items-center gap-2">
               <img src={logoImage} alt="Logo" className="h-8 w-8" />
               <span className="font-bold">PolyChemAI</span>
            </div>
            <button 
              onClick={() => setIsMobileMenuOpen(true)}
              className="p-2 bg-gray-100 dark:bg-slate-700 rounded-md text-gray-600 dark:text-gray-200"
            >
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h16" /></svg>
            </button>
        </div>

        <div className="p-4 md:p-8 min-h-full">
          {children}
        </div>
      </main>

      <ConfirmationModal
        isOpen={isLogoutModalOpen}
        title="Konfirmasi Keluar"
        message="Apakah Anda yakin ingin keluar?"
        onConfirm={handleConfirmLogout}
        onCancel={() => setIsLogoutModalOpen(false)}
        isLoading={isLoggingOut}
        isDanger={true}
      />
    </div>
  );
}

export default DashboardLayout;