import { Routes, Route } from "react-router-dom";
// 1. IMPORT TOASTER
import { Toaster } from "react-hot-toast";

import DashboardLayout from "./components/layout/DashboardLayout";
import HomePage from "./pages/HomePage";
import AboutPage from "./pages/AboutPage";
import ChemicalDetailPage from "./pages/ChemicalDetailPage";
import LoginPage from "./pages/LoginPage";
import RegisterPage from "./pages/RegisterPage";
import HistoryPage from "./pages/HistoryPage";
import LibraryPage from "./pages/LibraryPage";
import ForgotPasswordPage from "./pages/ForgotPasswordPage";
import WikiPage from "./pages/WikiPage";
import VerifyEmailNoticePage from "./pages/VerifyEmailNoticePage";
import { ThemeProvider } from "./context/ThemeContext";
// DITAMBAHKAN: guard rute untuk halaman yang wajib login
import ProtectedRoute from "./components/ProtectedRoute";

function App() {
  return (
    <ThemeProvider>
      {/* TOASTER */}
      <Toaster
        position="top-center"
        reverseOrder={false}
        toastOptions={{
          style: {
            background: "#1f2937", // Warna gelap (Dark Grey/Blue)
            color: "#fff", // Teks putih
            border: "1px solid #374151",
            padding: "16px",
            borderRadius: "8px",
          },
          success: {
            iconTheme: {
              primary: "#10b981", // Hijau bagus
              secondary: "#fff",
            },
          },
          error: {
            iconTheme: {
              primary: "#ef4444", // Merah bagus
              secondary: "#fff",
            },
          },
        }}
      />

      <Routes>
        {/* --- ZONE BEBAS (Tanpa Sidebar) --- */}
        <Route path="/login" element={<LoginPage />} />
        <Route path="/register" element={<RegisterPage />} />
        <Route path="/forgot-password" element={<ForgotPasswordPage />} />
        <Route path="/verify-email-notice" element={<VerifyEmailNoticePage />} />

        {/* --- ZONE DASHBOARD (Dengan Sidebar) --- */}
        <Route element={<DashboardLayout />}>
          {/* Halaman Utama — tetap publik, tamu boleh coba prediksi */}
          <Route path="/" element={<HomePage />} />
          {/* Halaman Detail — tetap publik, tamu boleh lihat hasil prediksi */}
          <Route path="/chemical/:id" element={<ChemicalDetailPage />} />
          {/* Halaman Output — tetap publik */}
          <Route path="/result" element={<ChemicalDetailPage />} />

          {/* Halaman History — DIUBAH: wajib login, redirect ke /login jika belum */}
          <Route
            path="/history"
            element={
              <ProtectedRoute>
                <HistoryPage />
              </ProtectedRoute>
            }
          />

          {/* Halaman Library — DIUBAH: wajib login, redirect ke /login jika belum */}
          <Route
            path="/library"
            element={
              <ProtectedRoute>
                <LibraryPage />
              </ProtectedRoute>
            }
          />

          {/* Halaman Wiki */}
          <Route path="/wiki" element={<WikiPage />} />
          {/* Halaman About */}
          <Route path="/about" element={<AboutPage />} />
        </Route>
      </Routes>
    </ThemeProvider>
  );
}

export default App;