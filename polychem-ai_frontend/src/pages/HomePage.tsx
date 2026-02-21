import { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { predictPolymer } from "../services/apiService";
import toast from "react-hot-toast";
import { auth } from "../lib/firebase";
import { addToHistory } from "../services/dbService";

function HomePage() {
  const [query, setQuery] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const navigate = useNavigate();

  // --- EFEK KETIK (TYPEWRITER) ---
  const textToType = "Plastic Discovery Agent";
  const [displayedText, setDisplayedText] = useState("");

  useEffect(() => {
    let index = 0;
    setDisplayedText("");

    const intervalId = setInterval(() => {
      setDisplayedText((prev) => textToType.slice(0, prev.length + 1));
      index++;
      if (index === textToType.length) clearInterval(intervalId);
    }, 100);

    return () => clearInterval(intervalId);
  }, []);

  // --- FUNGSI SEARCH ---
  const handleSearch = async () => {
    const cleanedQuery = query.trim();
    if (!cleanedQuery) {
      toast.error("Masukkan nama senyawa atau SMILES dulu!");
      return;
    }

    setIsLoading(true);

    try {
      const result = await predictPolymer(cleanedQuery);

      if (result) {
        toast.success("Prediksi berhasil ditemukan!");

        // 2. Simpan History (jika login)
        const user = auth.currentUser;
        if (user) {
          addToHistory(user.uid, cleanedQuery, result);
        }

        // 3. Pindah ke halaman Result
        navigate("/result", { state: { predictionData: result } });
      } else {
        toast.error("Input tidak valid / data tidak ditemukan.");
      }
    } catch (err) {
      console.error("Error detail:", err);
      toast.error("Gagal terhubung ke Backend.");
    } finally {
      setIsLoading(false);
    }
  };

  return (
    // Layout Utama: Flex Column, Centered, Full Height minus Header/Footer
    <div className="min-h-[70vh] md:min-h-[80vh] flex flex-col justify-center items-center pb-20 md:pb-40 max-w-5xl mx-auto px-4 md:px-8">
      {/* HEADER SECTION */}
      <div className="text-center mb-8 md:mb-10 w-full">
        {/* JUDUL dengan Efek Ketik */}
        <h1 className="text-3xl sm:text-4xl md:text-6xl font-extrabold text-main tracking-tight mb-4 md:mb-6 min-h-[3rem] md:min-h-[5rem] flex justify-center items-center">
          <span>{displayedText}</span>
          <span className="animate-pulse text-blue-500 ml-1">|</span>
        </h1>

        {/* SUBTITLE */}
        <p className="text-sm sm:text-base md:text-lg text-muted max-w-xs sm:max-w-2xl mx-auto leading-relaxed opacity-90">
          Masukkan struktur SMILES untuk memprediksi properti kimia dan jalur
          sintesis berbasis AI.
        </p>
      </div>

      {/* SEARCH BOX WRAPPER */}
      <div className="w-full max-w-3xl">
        <div className="bg-card p-1.5 md:p-2 rounded-xl md:rounded-2xl shadow-xl border border-border flex items-center transition-all duration-300 focus-within:shadow-2xl focus-within:border-blue-500 hover:shadow-2xl">
          {/* INPUT FIELD */}
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && handleSearch()}
            placeholder="Cari SMILES (Contoh: C=C)..."
            className="flex-1 px-3 md:px-5 py-3 md:py-4 text-main placeholder-gray-400 text-base md:text-lg bg-transparent outline-none w-full min-w-0"
            disabled={isLoading}
          />

          {/* TOMBOL SEARCH */}
          <button
            onClick={handleSearch}
            disabled={isLoading}
            className={`px-5 py-2.5 md:px-8 md:py-4 rounded-lg md:rounded-xl font-medium md:font-bold transition-all duration-300 shadow-md text-white text-sm md:text-base whitespace-nowrap
              ${
                isLoading
                  ? "bg-gray-400 cursor-not-allowed"
                  : "bg-slate-900 hover:bg-slate-800 dark:bg-blue-600 dark:hover:bg-blue-500 active:scale-95 md:hover:-translate-y-1"
              }`}
          >
            {isLoading ? "Running..." : "SEARCH"}
          </button>
        </div>
      </div>
    </div>
  );
}

export default HomePage;
