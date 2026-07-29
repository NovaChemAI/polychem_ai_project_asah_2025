import { useState, useEffect, useRef } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { createUserWithEmailAndPassword, updateProfile, sendEmailVerification } from 'firebase/auth';
import ReCAPTCHA from 'react-google-recaptcha';
import { auth } from '../lib/firebase';
import { syncUserProfile, verifyCaptcha } from '../services/apiService';
import { useTheme } from '../context/ThemeContext';

// ▼ IMPORT GAMBAR DARI ASSETS ▼
import registerBg from '../assets/RegisterImage.jpg'; 

function RegisterPage() {
  const navigate = useNavigate();
  const { theme } = useTheme();

  // State Input
  const [name, setName] = useState('');
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  
  // State UI & Validasi
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [isPasswordValid, setIsPasswordValid] = useState(false);

  // State Captcha
  const [captchaToken, setCaptchaToken] = useState<string | null>(null);
  const recaptchaRef = useRef<ReCAPTCHA>(null);

  // State fokus input password — untuk show/hide tooltip kriteria
  const [isPasswordFocused, setIsPasswordFocused] = useState(false);

  // State Kriteria Password
  const [passwordCriteria, setPasswordCriteria] = useState({
    length: false,
    upper: false,
    lower: false,
    number: false,
    special: false
  });

  // Validasi Real-time
  useEffect(() => {
    const criteria = {
      length: password.length >= 8,
      upper: /[A-Z]/.test(password),
      lower: /[a-z]/.test(password),
      number: /[0-9]/.test(password),
      special: /[!@#$%^&*(),.?":{}|<>_]/.test(password)
    };
    setPasswordCriteria(criteria);
    setIsPasswordValid(Object.values(criteria).every(Boolean));
  }, [password]);

  const handleRegister = async (e: React.FormEvent) => {
    e.preventDefault();

    if (!isPasswordValid) {
      setError("Mohon penuhi semua kriteria password.");
      return;
    }

    if (!captchaToken) {
      setError("Mohon selesaikan verifikasi CAPTCHA.");
      return;
    }

    setError('');
    setLoading(true);

    try {
      // 🔍 DEBUG SEMENTARA — hapus setelah masalah captcha selesai
      console.log('DEBUG captchaToken:', captchaToken);
      console.log('DEBUG panjang token:', captchaToken.length);

      // A. Verifikasi captcha ke backend dulu
      const captchaValid = await verifyCaptcha(captchaToken);
      if (!captchaValid) {
        setError('Verifikasi CAPTCHA gagal, silakan coba lagi.');
        recaptchaRef.current?.reset();
        setCaptchaToken(null);
        setLoading(false);
        return;
      }

      // B. Buat Akun (Firebase)
      const userCredential = await createUserWithEmailAndPassword(auth, email, password);
      const user = userCredential.user;

      // Update Profile
      await updateProfile(user, {
        displayName: name
      });

      // C. Kirim email verifikasi
      await sendEmailVerification(user);

      // D. Simpan ke Database
      await syncUserProfile(name);

      console.log("User registered & saved to DB:", user);

      // E. Arahkan ke halaman "cek email", bukan langsung ke home
      navigate('/verify-email-notice');

    } catch (error) {
      console.error(error);
      const err = error as { code?: string };

      if (err.code === 'auth/email-already-in-use') {
        setError('Email ini sudah terdaftar. Silakan login.');
      } else {
        setError('Gagal mendaftar. Silakan cek koneksi atau coba lagi.');
      }
      recaptchaRef.current?.reset();
      setCaptchaToken(null);
    } finally {
      setLoading(false);
    }
  };

  // Teks kriteria inline — bold, berubah hijau saat kriteria terpenuhi
  const RequirementText = ({ met, children }: { met: boolean, children: React.ReactNode }) => (
    <span className={`font-bold transition-colors duration-200 ${met ? 'text-green-600 dark:text-green-400' : 'text-main'}`}>
      {children}
    </span>
  );

  return (
    // Update: bg-background
    <div className="min-h-screen flex bg-background transition-colors duration-300">
      
      {/* KIRI: Ilustrasi */}
      <div className="hidden lg:block lg:w-1/2 relative">
        <img 
          src={registerBg} 
          alt="Lab Research" 
          className="absolute inset-0 h-full w-full object-cover"
        />
        {/* Update: Dark overlay for text contrast */}
        <div className="absolute inset-0 bg-blue-900/60 backdrop-blur-[2px] flex items-center justify-center">
          <div className="text-center p-10 text-white">
            <h2 className="text-4xl font-bold mb-4">Bergabung Bersama Kami</h2>
            <p className="text-lg text-blue-100">Mulai perjalanan penemuan senyawa barumu hari ini.</p>
          </div>
        </div>
      </div>

      {/* KANAN: Form Register */}
      <div className="w-full lg:w-1/2 flex flex-col justify-center px-8 py-12 sm:px-12 lg:px-16 xl:px-24">
        <div className="sm:mx-auto sm:w-full sm:max-w-md mb-8">
          {/* Logo tetap ambil dari public */}
          <img className="h-12 w-auto" src="/images/logo.png" alt="NovaChem Logo" />
          {/* Update: text-main */}
          <h2 className="mt-6 text-3xl font-extrabold text-main">
            Buat Akun Baru
          </h2>
          {/* Update: text-muted */}
          <p className="mt-2 text-sm text-muted">
            Sudah punya akun?{' '}
            {/* Update: Link color dark mode */}
            <Link to="/login" className="font-medium text-blue-600 dark:text-blue-400 hover:text-blue-500 hover:underline">
              Login di sini
            </Link>
          </p>
        </div>

        <div className="sm:mx-auto sm:w-full sm:max-w-md">
          
          {error && (
            // Update: Dark mode error colors
            <div className="mb-4 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 text-red-600 dark:text-red-400 px-4 py-3 rounded-lg text-sm">
              {error}
            </div>
          )}

          <form className="space-y-5" onSubmit={handleRegister}>
            <div>
              {/* Update: text-main */}
              <label className="block text-sm font-medium text-main">Nama Lengkap</label>
              <input 
                type="text" 
                required 
                // Update: bg-input-bg, border-border, text-main
                className="mt-1 block w-full px-3 py-3 border border-border rounded-lg shadow-sm focus:ring-blue-500 focus:border-blue-500 sm:text-sm bg-input-bg text-main" 
                value={name} 
                onChange={(e) => setName(e.target.value)}
              />
            </div>
            
            <div>
              <label className="block text-sm font-medium text-main">Email address</label>
              <input 
                type="email" 
                required 
                className="mt-1 block w-full px-3 py-3 border border-border rounded-lg shadow-sm focus:ring-blue-500 focus:border-blue-500 sm:text-sm bg-input-bg text-main" 
                value={email} 
                onChange={(e) => setEmail(e.target.value)}
              />
            </div>
            
            <div>
              <label className="block text-sm font-medium text-main">Password</label>
              <div className="relative mt-1">
                <input 
                  type="password" 
                  required 
                  // Update: conditional border colors for dark mode
                  className={`block w-full px-3 py-3 pr-10 border rounded-lg shadow-sm sm:text-sm bg-input-bg text-main 
                    ${isPasswordValid 
                      ? 'border-green-500 focus:ring-green-500 dark:border-green-400' 
                      : 'border-border focus:ring-blue-500'}`} 
                  value={password} 
                  onChange={(e) => setPassword(e.target.value)}
                  onFocus={() => setIsPasswordFocused(true)}
                  onBlur={() => setIsPasswordFocused(false)}
                />
                {isPasswordValid && (
                  <svg
                    className="absolute right-3 top-1/2 -translate-y-1/2 w-5 h-5 text-green-600 dark:text-green-400"
                    fill="none" stroke="currentColor" strokeWidth={3} viewBox="0 0 24 24"
                  >
                    <path strokeLinecap="round" strokeLinejoin="round" d="M4.5 12.75l6 6 9-13.5" />
                  </svg>
                )}
              </div>

              {/* Tooltip Kriteria Password — hanya muncul saat input di-fokus */}
              <div
                className={`overflow-hidden transition-all duration-200 ease-out ${
                  isPasswordFocused ? 'max-h-40 opacity-100 pt-4' : 'max-h-0 opacity-0 pt-0'
                }`}
              >
                <div className="relative">
                  {/* Anak panah menunjuk ke input */}
                  <div className="absolute -top-[7px] left-6 w-3 h-3 bg-gray-50 dark:bg-slate-800/60 border-l border-t border-border rotate-45" />
                  <div className="p-3.5 rounded-lg border border-border bg-gray-50 dark:bg-slate-800/60 text-sm leading-relaxed text-main">
                    Password harus minimal <RequirementText met={passwordCriteria.length}>8 karakter</RequirementText> dan mengandung{' '}
                    <RequirementText met={passwordCriteria.number}>1 angka</RequirementText>,{' '}
                    <RequirementText met={passwordCriteria.upper}>1 huruf besar</RequirementText>,{' '}
                    <RequirementText met={passwordCriteria.lower}>1 huruf kecil</RequirementText>, dan{' '}
                    <RequirementText met={passwordCriteria.special}>1 simbol</RequirementText>.
                  </div>
                </div>
              </div>
            </div>

            {/* Widget CAPTCHA */}
            <div className="rounded-xl border border-border bg-gray-50 dark:bg-slate-800/60 overflow-hidden">
              <div className="flex items-center justify-between px-4 pt-3.5">
                <span className="text-[11px] font-semibold uppercase tracking-wide text-gray-400 dark:text-slate-500">
                  Verifikasi Keamanan
                </span>
                {captchaToken && (
                  <span className="inline-flex items-center gap-1 text-[11px] font-semibold text-green-600 dark:text-green-400">
                    <svg className="w-3.5 h-3.5" fill="none" stroke="currentColor" strokeWidth={3} viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" d="M4.5 12.75l6 6 9-13.5" />
                    </svg>
                    Terverifikasi
                  </span>
                )}
              </div>
              <div className="flex justify-center p-3.5 pt-2.5">
                <div className="rounded-lg overflow-hidden shadow-sm">
                  <ReCAPTCHA
                    key={theme}
                    ref={recaptchaRef}
                    sitekey={import.meta.env.VITE_RECAPTCHA_SITE_KEY}
                    theme={theme}
                    onChange={(token) => setCaptchaToken(token)}
                    onExpired={() => setCaptchaToken(null)}
                  />
                </div>
              </div>
            </div>
            
            <div className="pt-2">
                <button 
                  type="submit" 
                  disabled={loading || !isPasswordValid || !captchaToken} 
                  // Update: Button colors (Slate/Blue)
                  className={`w-full flex justify-center py-3 px-4 border border-transparent rounded-lg shadow-sm text-sm font-bold text-white transition-colors
                    ${loading || !isPasswordValid || !captchaToken
                      ? 'bg-gray-400 dark:bg-slate-700 cursor-not-allowed' 
                      : 'bg-slate-900 dark:bg-blue-600 hover:bg-slate-800 dark:hover:bg-blue-700'}`}
                >
                  {loading ? 'Mendaftarkan...' : 'Daftar Sekarang'}
                </button>
            </div>
          </form>
        </div>
      </div>
    </div>
  );
}

export default RegisterPage;