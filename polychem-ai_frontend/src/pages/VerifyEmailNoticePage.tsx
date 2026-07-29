import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { sendEmailVerification, signOut } from 'firebase/auth';
import toast from 'react-hot-toast';
import { auth } from '../lib/firebase';

const RESEND_COOLDOWN_SECONDS = 60;

function VerifyEmailNoticePage() {
  const navigate = useNavigate();
  const [sending, setSending] = useState(false);
  const [cooldown, setCooldown] = useState(0);
  const [checking, setChecking] = useState(false);

  const email = auth.currentUser?.email ?? '';

  // Hitung mundur cooldown tombol kirim ulang
  useEffect(() => {
    if (cooldown <= 0) return;
    const timer = setInterval(() => {
      setCooldown((prev) => (prev <= 1 ? 0 : prev - 1));
    }, 1000);
    return () => clearInterval(timer);
  }, [cooldown]);

  const handleResend = async () => {
    if (!auth.currentUser) {
      toast.error('Sesi tidak ditemukan. Silakan login ulang.');
      navigate('/login');
      return;
    }

    setSending(true);
    try {
      await sendEmailVerification(auth.currentUser);
      toast.success('Email verifikasi telah dikirim ulang!');
      setCooldown(RESEND_COOLDOWN_SECONDS);
    } catch (error) {
      console.error(error);
      toast.error('Gagal mengirim ulang email. Coba lagi beberapa saat.');
    } finally {
      setSending(false);
    }
  };

  // Cek apakah user sudah verifikasi (reload data user dari Firebase)
  const handleCheckVerified = async () => {
    if (!auth.currentUser) {
      navigate('/login');
      return;
    }

    setChecking(true);
    try {
      await auth.currentUser.reload();
      if (auth.currentUser.emailVerified) {
        toast.success('Email berhasil diverifikasi!');
        navigate('/');
      } else {
        toast.error('Email belum diverifikasi. Silakan cek inbox kamu.');
      }
    } catch (error) {
      console.error(error);
      toast.error('Gagal memeriksa status verifikasi.');
    } finally {
      setChecking(false);
    }
  };

  const handleLogout = async () => {
    await signOut(auth);
    navigate('/login');
  };

  return (
    <div className="min-h-screen flex items-center justify-center bg-background px-4 transition-colors duration-300">
      <div className="w-full max-w-md">
        <div className="bg-white dark:bg-slate-800/60 border border-border rounded-2xl shadow-sm p-8 text-center">
          {/* Ikon amplop */}
          <div className="mx-auto mb-5 w-16 h-16 rounded-full bg-blue-50 dark:bg-blue-900/30 flex items-center justify-center">
            <svg
              className="w-8 h-8 text-blue-600 dark:text-blue-400"
              fill="none" stroke="currentColor" strokeWidth={1.75} viewBox="0 0 24 24"
            >
              <path strokeLinecap="round" strokeLinejoin="round" d="M21.75 6.75v10.5a2.25 2.25 0 01-2.25 2.25h-15a2.25 2.25 0 01-2.25-2.25V6.75m19.5 0A2.25 2.25 0 0019.5 4.5h-15a2.25 2.25 0 00-2.25 2.25m19.5 0v.243a2.25 2.25 0 01-1.07 1.916l-7.5 4.615a2.25 2.25 0 01-2.36 0L3.32 8.91a2.25 2.25 0 01-1.07-1.916V6.75" />
            </svg>
          </div>

          <h1 className="text-2xl font-extrabold text-main mb-2">
            Cek Email Kamu
          </h1>
          <p className="text-sm text-muted leading-relaxed">
            Kami sudah mengirim link verifikasi ke{' '}
            <span className="font-semibold text-main">{email || 'email kamu'}</span>.
            Klik link tersebut untuk mengaktifkan akun.
          </p>

          <div className="mt-6 space-y-3">
            <button
              onClick={handleCheckVerified}
              disabled={checking}
              className={`w-full flex justify-center py-3 px-4 rounded-lg shadow-sm text-sm font-bold text-white transition-colors
                ${checking
                  ? 'bg-gray-400 dark:bg-slate-700 cursor-not-allowed'
                  : 'bg-slate-900 dark:bg-blue-600 hover:bg-slate-800 dark:hover:bg-blue-700'}`}
            >
              {checking ? 'Memeriksa...' : 'Saya Sudah Verifikasi'}
            </button>

            <button
              onClick={handleResend}
              disabled={sending || cooldown > 0}
              className={`w-full flex justify-center py-3 px-4 rounded-lg border text-sm font-semibold transition-colors
                ${sending || cooldown > 0
                  ? 'border-border text-gray-400 dark:text-slate-500 cursor-not-allowed'
                  : 'border-border text-main hover:bg-gray-50 dark:hover:bg-slate-700/50'}`}
            >
              {cooldown > 0
                ? `Kirim ulang dalam ${cooldown}s`
                : sending
                ? 'Mengirim...'
                : 'Kirim Ulang Email Verifikasi'}
            </button>
          </div>

          <button
            onClick={handleLogout}
            className="mt-5 text-xs text-muted hover:text-main hover:underline transition-colors"
          >
            Bukan kamu? Keluar dan login dengan akun lain
          </button>
        </div>
      </div>
    </div>
  );
}

export default VerifyEmailNoticePage;