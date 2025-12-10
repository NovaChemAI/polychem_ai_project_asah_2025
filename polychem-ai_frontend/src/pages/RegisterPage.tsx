import { useState, useEffect } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { createUserWithEmailAndPassword, updateProfile } from 'firebase/auth';
import { doc, setDoc } from 'firebase/firestore'; 
import { auth, db } from '../lib/firebase'; 

// ▼ IMPORT GAMBAR DARI ASSETS ▼
import registerBg from '../assets/RegisterImage.jpg'; 

function RegisterPage() {
  const navigate = useNavigate();

  // State Input
  const [name, setName] = useState('');
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  
  // State UI & Validasi
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const [isPasswordValid, setIsPasswordValid] = useState(false);

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

    setError('');
    setLoading(true);

    try {
      // A. Buat Akun
      const userCredential = await createUserWithEmailAndPassword(auth, email, password);
      const user = userCredential.user;

      // B. Update Profile
      await updateProfile(user, {
        displayName: name
      });

      // C. Simpan ke Database
      await setDoc(doc(db, "users", user.uid), {
        uid: user.uid,
        name: name,
        email: email,
        role: "user",
        createdAt: new Date(),
        searchHistory: [] 
      });

      console.log("User registered & saved to DB:", user);
      navigate('/');
      
    } catch (error) {
      console.error(error);
      const err = error as { code: string };

      if (err.code === 'auth/email-already-in-use') {
        setError('Email ini sudah terdaftar. Silakan login.');
      } else {
        setError('Gagal mendaftar. Silakan cek koneksi atau coba lagi.');
      }
    } finally {
      setLoading(false);
    }
  };

  // Komponen Checklist
  const RequirementItem = ({ met, text }: { met: boolean, text: string }) => (
    <div className={`flex items-center text-xs mt-1 ${met ? 'text-green-600' : 'text-gray-400'}`}>
      <span className={`mr-2 flex-shrink-0 w-4 h-4 flex items-center justify-center rounded-full border ${met ? 'border-green-600 bg-green-100' : 'border-gray-300'}`}>
        {met && <svg className="w-2.5 h-2.5" fill="currentColor" viewBox="0 0 20 20"><path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd"/></svg>}
      </span>
      {text}
    </div>
  );

  return (
    <div className="min-h-screen flex bg-white">
      {/* KIRI: Ilustrasi */}
      <div className="hidden lg:block lg:w-1/2 relative">
        {}
        <img 
          src={registerBg} 
          alt="Lab Research" 
          className="absolute inset-0 h-full w-full object-cover"
        />
        <div className="absolute inset-0 bg-blue-900/40 backdrop-blur-[2px] flex items-center justify-center">
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
          <h2 className="mt-6 text-3xl font-extrabold text-gray-900">
            Buat Akun Baru
          </h2>
          <p className="mt-2 text-sm text-gray-600">
            Sudah punya akun?{' '}
            <Link to="/login" className="font-medium text-blue-600 hover:text-blue-500">
              Login di sini
            </Link>
          </p>
        </div>

        <div className="sm:mx-auto sm:w-full sm:max-w-md">
          
          {error && <div className="mb-4 bg-red-50 border border-red-200 text-red-600 px-4 py-3 rounded-lg text-sm">{error}</div>}

          <form className="space-y-5" onSubmit={handleRegister}>
            <div>
              <label className="block text-sm font-medium text-gray-700">Nama Lengkap</label>
              <input type="text" required className="mt-1 block w-full px-3 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-blue-500 focus:border-blue-500 sm:text-sm" value={name} onChange={(e) => setName(e.target.value)}/>
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700">Email address</label>
              <input type="email" required className="mt-1 block w-full px-3 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-blue-500 focus:border-blue-500 sm:text-sm" value={email} onChange={(e) => setEmail(e.target.value)}/>
            </div>
            <div>
              <label className="block text-sm font-medium text-gray-700">Password</label>
              <input type="password" required className={`mt-1 block w-full px-3 py-3 border rounded-lg shadow-sm sm:text-sm ${isPasswordValid ? 'border-green-500 focus:ring-green-500' : 'border-gray-300 focus:ring-blue-500'}`} value={password} onChange={(e) => setPassword(e.target.value)}/>
              <div className="mt-3 bg-gray-50 p-3 rounded-lg border border-gray-100 grid grid-cols-1 sm:grid-cols-2 gap-1">
                <RequirementItem met={passwordCriteria.length} text="Minimal 8 karakter" />
                <RequirementItem met={passwordCriteria.upper} text="Huruf Besar (A-Z)" />
                <RequirementItem met={passwordCriteria.lower} text="Huruf Kecil (a-z)" />
                <RequirementItem met={passwordCriteria.number} text="Angka (0-9)" />
                <RequirementItem met={passwordCriteria.special} text="Simbol (!@#_)" />
              </div>
            </div>
            <div className="pt-2">
                <button type="submit" disabled={loading || !isPasswordValid} className={`w-full flex justify-center py-3 px-4 border border-transparent rounded-lg shadow-sm text-sm font-medium text-white ${loading || !isPasswordValid ? 'bg-gray-400 cursor-not-allowed' : 'bg-blue-600 hover:bg-blue-700'}`}>{loading ? 'Mendaftarkan...' : 'Daftar Sekarang'}</button>
            </div>
          </form>
        </div>
      </div>
    </div>
  );
}

export default RegisterPage;