import { Routes, Route } from 'react-router-dom';
import DashboardLayout from './components/layout/DashboardLayout';
import HomePage from './pages/HomePage';
import ExplorePage from './pages/ExplorePage';
import AboutPage from './pages/AboutPage';
import ChemicalDetailPage from './pages/ChemicalDetailPage';
import LoginPage from './pages/LoginPage';
import RegisterPage from './pages/RegisterPage';
import HistoryPage from './pages/HistoryPage';
import LibraryPage from './pages/LibraryPage';
import ForgotPasswordPage from './pages/ForgotPasswordPage'; // 1. Pastikan Import Ini Ada
import WikiPage from './pages/WikiPage';

function App() {
  return (
    <Routes>
      {/* --- ZONE BEBAS (Tanpa Sidebar) --- */}
      {/* Rute ini harus ditaruh PALING ATAS agar dicek duluan */}
      
      <Route path="/login" element={<LoginPage />} />
      <Route path="/register" element={<RegisterPage />} />
      
      {/* ▼▼▼ POSISI RUTE INI HARUS DI SINI ▼▼▼ */}
      <Route path="/forgot-password" element={<ForgotPasswordPage />} />
      {/* ▲▲▲ JANGAN TARUH DI BAWAH ATAU DI DALAM DASHBOARD ▲▲▲ */}


      {/* --- ZONE DASHBOARD (Dengan Sidebar) --- */}
      {/* path="*" artinya: "Kalau tidak cocok dengan yang di atas, masuk ke sini" */}
      <Route path="*" element={
        <DashboardLayout>
          <Routes>
            <Route path="/" element={<HomePage />} />
            <Route path="/explore" element={<ExplorePage />} />
            <Route path="/chemical/:id" element={<ChemicalDetailPage />} />
            <Route path="/history" element={<HistoryPage />} />
            <Route path="/library" element={<LibraryPage />} />
            <Route path="/wiki" element={<WikiPage />} />
            <Route path="/about" element={<AboutPage />} />
          </Routes>
        </DashboardLayout>
      } />
    </Routes>
  )
}

export default App;