import { useState } from 'react';
import { Link } from 'react-router-dom';
import { sendPasswordResetEmail } from 'firebase/auth'; 
import { auth } from '../lib/firebase';
import logoImage from '/images/logo.png';
import loginBg from '../assets/LoginImage.jpg'; 

function ForgotPasswordPage() {
  const [email, setEmail] = useState('');
  const [message, setMessage] = useState(''); 
  const [error, setError] = useState('');    
  const [loading, setLoading] = useState(false);

  const handleResetPassword = async (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!email) {
      setError("Mohon masukkan email Anda.");
      return;
    }

    setError('');
    setMessage('');
    setLoading(true);

    try {
      // --- LOGIC FIREBASE ---
      await sendPasswordResetEmail(auth, email);
      
      // Jika berhasil
      setMessage('Link reset password telah dikirim ke email Anda. Silakan cek Inbox atau Spam.');
      setEmail(''); // Kosongkan form
    } catch (error) {
      console.error(error);
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const err = error as { code: string };
      
      if (err.code === 'auth/user-not-found') {
        setError('Email ini tidak terdaftar di sistem kami.');
      } else if (err.code === 'auth/invalid-email') {
        setError('Format email tidak valid.');
      } else {
        setError('Gagal mengirim email. Coba lagi nanti.');
      }
    } finally {
      setLoading(false);
    }
  };

  return (
    // Update: bg-white -> bg-background
    <div className="min-h-screen flex bg-background transition-colors duration-300">
      
      {/* BAGIAN KIRI (GAMBAR) */}
      <div className="hidden lg:block lg:w-1/2 relative">
        <img 
          src={loginBg} 
          alt="Chemistry Lab" 
          className="absolute inset-0 h-full w-full object-cover"
        />
        {/* Overlay tetap gelap agar teks putih terbaca */}
        <div className="absolute inset-0 bg-slate-900/60 backdrop-blur-[2px] flex items-center justify-center">
          <div className="text-center p-10 text-white">
            <h2 className="text-4xl font-bold mb-4">Recovery Akses</h2>
            <p className="text-lg text-blue-100">Kembalikan akses akun Anda untuk melanjutkan riset.</p>
          </div>
        </div>
      </div>

      {/* BAGIAN KANAN (FORM) */}
      <div className="w-full lg:w-1/2 flex flex-col justify-center px-8 py-12 sm:px-12 lg:px-16 xl:px-24">
        
        <div className="sm:mx-auto sm:w-full sm:max-w-md mb-8">
          <img className="h-12 w-auto" src={logoImage} alt="NovaChem Logo" />
          {/* Update: text-gray-900 -> text-main */}
          <h2 className="mt-6 text-3xl font-extrabold text-main">
            Lupa Password?
          </h2>
          {/* Update: text-gray-600 -> text-muted */}
          <p className="mt-2 text-sm text-muted">
            Jangan khawatir. Masukkan email Anda dan kami akan mengirimkan petunjuk reset password.
          </p>
        </div>

        <div className="sm:mx-auto sm:w-full sm:max-w-md">
          
          {/* Notifikasi Sukses (Hijau) */}
          {message && (
            // Update: Dark mode variant colors
            <div className="mb-4 bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800 text-green-700 dark:text-green-400 px-4 py-3 rounded-lg text-sm flex items-start">
              <svg className="h-5 w-5 mr-2 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M9 12l2 2 4-4m6 2a9 9 0 11-18 0 9 9 0 0118 0z" /></svg>
              {message}
            </div>
          )}

          {/* Notifikasi Error (Merah) */}
          {error && (
            // Update: Dark mode variant colors
            <div className="mb-4 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 text-red-600 dark:text-red-400 px-4 py-3 rounded-lg text-sm">
              {error}
            </div>
          )}

          <form className="space-y-6" onSubmit={handleResetPassword}>
            <div>
              {/* Update: text-main */}
              <label htmlFor="email" className="block text-sm font-medium text-main">Email address</label>
              <div className="mt-1">
                <input
                  id="email"
                  type="email"
                  required
                  // Update: bg-input-bg, border-border, text-main
                  className="appearance-none block w-full px-3 py-3 border border-border rounded-lg shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 sm:text-sm bg-input-bg text-main"
                  placeholder="contoh@email.com"
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                />
              </div>
            </div>

            <div>
              <button
                type="submit"
                disabled={loading}
                // Update: Button color consistent with DetailPage (Slate in light, Blue in dark)
                className="w-full flex justify-center py-3 px-4 border border-transparent rounded-lg shadow-sm text-sm font-bold text-white bg-slate-900 dark:bg-blue-600 hover:bg-slate-800 dark:hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 disabled:opacity-70 disabled:cursor-not-allowed transition-colors"
              >
                {loading ? 'Mengirim...' : 'Kirim Link Reset'}
              </button>
            </div>
          </form>

          <div className="mt-6 text-center">
            {/* Update: text-blue-600 dark:text-blue-400 */}
            <Link to="/login" className="font-medium text-blue-600 dark:text-blue-400 hover:text-blue-500 dark:hover:text-blue-300 flex items-center justify-center transition-colors">
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M10 19l-7-7m0 0l7-7m-7 7h18" /></svg>
              Kembali ke Login
            </Link>
          </div>

        </div>
      </div>
    </div>
  );
}

export default ForgotPasswordPage;