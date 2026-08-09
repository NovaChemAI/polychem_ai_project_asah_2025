import { useEffect, useState } from "react";
import { Navigate, useLocation } from "react-router-dom";
import { onAuthStateChanged } from "firebase/auth";
import { auth } from "../lib/firebase";

function ProtectedRoute({ children }: { children: React.ReactNode }) {
  const location = useLocation();
  const [checking, setChecking] = useState(true);
  const [isLoggedIn, setIsLoggedIn] = useState(false);

  useEffect(() => {
    const unsubscribe = onAuthStateChanged(auth, (user) => {
      setIsLoggedIn(!!user);
      setChecking(false);
    });
    return () => unsubscribe();
  }, []);

  if (checking) {
    // Sebentar menunggu Firebase memuat status auth, hindari "flash" redirect
    // yang salah sebelum status login benar-benar diketahui.
    return null; // atau ganti dengan komponen loading spinner Anda
  }

  if (!isLoggedIn) {
    // state: from → supaya nanti bisa diarahkan balik ke halaman asal setelah login (opsional, lihat poin 4)
    return <Navigate to="/login" state={{ from: location.pathname }} replace />;
  }

  return <>{children}</>;
}

export default ProtectedRoute;