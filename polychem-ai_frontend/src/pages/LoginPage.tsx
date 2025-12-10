import { useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { signInWithEmailAndPassword, signInWithPopup } from 'firebase/auth';
// Import Firestore untuk cek/simpan user Google
import { doc, getDoc, setDoc } from 'firebase/firestore'; 
import { auth, googleProvider, db } from '../lib/firebase';

// Import Gambar Background
import loginBg from '../assets/LoginImage.jpg'; 

function LoginPage() {
  const navigate = useNavigate();
  
  // State
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);

  // --- LOGIC 1: LOGIN EMAIL ---
  const handleLogin = async (e: React.FormEvent) => {
    e.preventDefault();
    setError('');
    setLoading(true);

    try {
      await signInWithEmailAndPassword(auth, email, password);
      navigate('/');
    } catch (error) {
      console.error(error);
      const err = error as { code: string }; 
      if (err.code === 'auth/invalid-credential') {
        setError('Email atau password salah.');
      } else {
        setError('Gagal login. Coba lagi nanti.');
      }
    } finally {
      setLoading(false);
    }
  };

  // --- LOGIC 2: LOGIN GOOGLE (Dengan Simpan Database) ---
  const handleGoogleLogin = async () => {
    setError(''); 
    try {
      // 1. Buka Popup Google
      const result = await signInWithPopup(auth, googleProvider);
      const user = result.user;

      // 2. Cek apakah user ini sudah ada di Database kita?
      const userDocRef = doc(db, "users", user.uid);
      const userDocSnap = await getDoc(userDocRef);

      // 3. Jika BELUM ada (User Baru), simpan datanya!
      if (!userDocSnap.exists()) {
        await setDoc(userDocRef, {
          uid: user.uid,
          name: user.displayName,
          email: user.email,
          role: "user",
          createdAt: new Date(),
          searchHistory: [],
          photoURL: user.photoURL
        });
        console.log("User Google baru disimpan ke DB");
      }

      // 4. Masuk Dashboard
      navigate('/');

    } catch (error) {
      console.error("Google Login Error:", error);
      setError('Gagal login dengan Google. Pastikan koneksi aman.');
    }
  };

  return (
    <div className="min-h-screen flex bg-white">
      
      {/* BAGIAN KIRI (GAMBAR) */}
      <div className="hidden lg:block lg:w-1/2 relative">
        <img 
          src={loginBg} 
          alt="Chemistry Lab" 
          className="absolute inset-0 h-full w-full object-cover"
        />
        <div className="absolute inset-0 bg-blue-900/40 backdrop-blur-[2px] flex items-center justify-center">
          <div className="text-center p-10 text-white">
            <h2 className="text-4xl font-bold mb-4">PolyChemAI</h2>
            <p className="text-lg text-blue-100">Mempercepat Penemuan Kimia dengan Kecerdasan Buatan.</p>
          </div>
        </div>
      </div>

      {/* BAGIAN KANAN (FORM) */}
      <div className="w-full lg:w-1/2 flex flex-col justify-center px-8 py-12 sm:px-12 lg:px-16 xl:px-24">
        
        <div className="sm:mx-auto sm:w-full sm:max-w-md mb-8">
          {/* Logo mengambil dari folder public */}
          <img className="h-12 w-auto" src="/images/logo.png" alt="PolyChem Logo" />
          <h2 className="mt-6 text-3xl font-extrabold text-gray-900">
            Selamat Datang Kembali
          </h2>
          <p className="mt-2 text-sm text-gray-600">
            Belum punya akun?{' '}
            <Link to="/register" className="font-medium text-blue-600 hover:text-blue-500">
              Daftar gratis di sini
            </Link>
          </p>
        </div>

        <div className="sm:mx-auto sm:w-full sm:max-w-md">
          
          {/* Pesan Error */}
          {error && (
            <div className="mb-4 bg-red-50 border border-red-200 text-red-600 px-4 py-3 rounded-lg text-sm">
              {error}
            </div>
          )}

          <form className="space-y-6" onSubmit={handleLogin}>
            <div>
              <label htmlFor="email" className="block text-sm font-medium text-gray-700">Email address</label>
              <div className="mt-1">
                <input
                  id="email"
                  type="email"
                  required
                  className="appearance-none block w-full px-3 py-3 border border-gray-300 rounded-lg shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 sm:text-sm"
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                />
              </div>
            </div>

            <div>
              {/* Label Password dengan Link Lupa Password */}
              <div className="flex items-center justify-between">
                <label htmlFor="password" className="block text-sm font-medium text-gray-700">
                  Password
                </label>
                <div className="text-sm">
                  <Link to="/forgot-password" className="font-medium text-blue-600 hover:text-blue-500">
                    Lupa password?
                  </Link>
                </div>
              </div>
              <div className="mt-1">
                <input
                  id="password"
                  type="password"
                  required
                  className="appearance-none block w-full px-3 py-3 border border-gray-300 rounded-lg shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 sm:text-sm"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                />
              </div>
            </div>

            <div>
              <button
                type="submit"
                disabled={loading}
                className="w-full flex justify-center py-3 px-4 border border-transparent rounded-lg shadow-sm text-sm font-medium text-white bg-blue-600 hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 disabled:bg-blue-300"
              >
                {loading ? 'Memproses...' : 'Sign in'}
              </button>
            </div>
          </form>

          <div className="mt-6">
            <div className="relative">
              <div className="absolute inset-0 flex items-center">
                <div className="w-full border-t border-gray-300" />
              </div>
              <div className="relative flex justify-center text-sm">
                <span className="px-2 bg-white text-gray-500">Atau masuk dengan</span>
              </div>
            </div>

            <div className="mt-6 grid grid-cols-2 gap-3">
              {/* Tombol Google */}
              <button
                onClick={handleGoogleLogin}
                type="button"
                className="w-full inline-flex justify-center items-center py-3 px-4 border border-gray-300 rounded-lg shadow-sm bg-white text-sm font-medium text-gray-500 hover:bg-gray-50 transition-colors"
              >
                <svg className="h-5 w-5 mr-2" viewBox="0 0 24 24"><path d="M22.56 12.25c0-.78-.07-1.53-.2-2.25H12v4.26h5.92c-.26 1.37-1.04 2.53-2.21 3.31v2.77h3.57c2.08-1.92 3.28-4.74 3.28-8.09z" fill="#4285F4"/><path d="M12 23c2.97 0 5.46-.98 7.28-2.66l-3.57-2.77c-.98.66-2.23 1.06-3.71 1.06-2.86 0-5.29-1.93-6.16-4.53H2.18v2.84C3.99 20.53 7.7 23 12 23z" fill="#34A853"/><path d="M5.84 14.09c-.22-.66-.35-1.36-.35-2.09s.13-1.43.35-2.09V7.07H2.18C1.43 8.55 1 10.22 1 12s.43 3.45 1.18 4.93l2.85-2.22.81-.62z" fill="#FBBC05"/><path d="M12 5.38c1.62 0 3.06.56 4.21 1.64l3.15-3.15C17.45 2.09 14.97 1 12 1 7.7 1 3.99 3.47 2.18 7.07l3.66 2.84c.87-2.6 3.3-4.53 6.16-4.53z" fill="#EA4335"/></svg>
                Google
              </button>

              {/* Tombol Tamu */}
              <button
                onClick={() => navigate('/')}
                type="button"
                className="w-full inline-flex justify-center items-center py-3 px-4 border border-gray-300 rounded-lg shadow-sm bg-white text-sm font-medium text-gray-500 hover:bg-gray-50 transition-colors"
              >
                <span className="mr-2">👤</span>
                Tamu
              </button>
            </div>
          </div>

        </div>
      </div>
    </div>
  );
}

export default LoginPage;