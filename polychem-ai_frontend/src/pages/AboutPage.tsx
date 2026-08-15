import { useState, useEffect } from 'react';
import reactLogo from '../assets/TechLogo/LogoReact.png';
import firebaseLogo from '../assets/TechLogo/LogoFirebase.png';
import pythonLogo from '../assets/TechLogo/LogoPhyton.png'; 
import tailwindLogo from '../assets/TechLogo/LogoTailwind.png';
import colabLogo from '../assets/TechLogo/LogoColab.png';
import vscodeLogo from '../assets/TechLogo/LogoVsCode.png';

// ▼ IMPORT GAMBAR ILUSTRASI TIM ▼
import userMale from '../assets/DevImages/ilustration user.jpg';
import userFemale from '../assets/DevImages/ilustration userw.png';

// 1. DEFINISIKAN TIPE DATA MEMBER
interface TeamMember {
  name: string;
  role: string;
  image: string;
  color: string;
}

function AboutPage() {
  const [isLoading, setIsLoading] = useState(true);

  useEffect(() => {
    const timer = setTimeout(() => setIsLoading(false), 1500); 
    return () => clearTimeout(timer);
  }, []);

  const technologies = [
    { name: "React + Vite", role: "Frontend Modern", image: reactLogo },
    { name: "Firebase", role: "Auth & Database", image: firebaseLogo },
    { name: "Python + RDKit", role: "Validasi & Fitur Kimia", image: pythonLogo },
    { name: "Tailwind CSS", role: "Styling System", image: tailwindLogo },
    { name: "Google Colab", role: "Eksperimen & Riset", image: colabLogo },
    { name: "VS Code", role: "Code Editor", image: vscodeLogo },
  ];

  const mlTeam: TeamMember[] = [
    { name: "Tiara Diansyah Putri", role: "Machine Learning", image: userFemale, color: "green" },
    { name: "Loista Amanda Noviar", role: "Machine Learning", image: userFemale, color: "green" }
  ];

  // --- EDIT BAGIAN INI ---
  const webTeam: TeamMember[] = [
    // Nama M Ihsan Sadik sudah dihapus dari sini
    { name: "Dikky Juliyanto", role: "Front-End & Back-End", image: userMale, color: "blue" },
    { name: "Zakiul Fata", role: "Front-End & Back-End", image: userMale, color: "blue" }
  ];

  // --- BARU: DATA PEMBIMBING ---
  const advisorTeam: TeamMember[] = [
    { name: "Nazaruddin Ahmad, M.T.", role: "Pembimbing 1", image: userMale, color: "blue" },
    { name: "Khairan AR, M.Kom", role: "Pembimbing 2", image: userMale, color: "blue" }
  ];

  const renderMemberCard = (member: TeamMember, index: number) => (
    <div key={index} className="bg-card p-6 rounded-2xl border border-border hover:shadow-lg transition-all group w-full sm:w-64 flex flex-col items-center text-center">
      <div className={`w-24 h-24 rounded-full mb-4 p-1 shadow-sm ring-4 ring-gray-50 dark:ring-slate-800 overflow-hidden bg-card
        ${member.color === 'green' ? 'group-hover:ring-green-100 dark:group-hover:ring-green-900' : 'group-hover:ring-blue-100 dark:group-hover:ring-blue-900'} transition-all`}>
        <img 
          src={member.image} 
          alt={member.name} 
          className="w-full h-full object-cover rounded-full"
        />
      </div>
      {/*▼ EDIT DI BAGIAN INI ▼*/}
      <h3 className="font-bold text-lg text-main mb-1 whitespace-nowrap truncate">{member.name}</h3>
      {/*▲ EDIT DI BAGIAN INI ▲*/}
      <span className={`px-3 py-1 rounded-full text-xs font-medium 
        ${member.color === 'green' 
          ? 'bg-green-50 text-green-700 dark:bg-green-900/30 dark:text-green-400' 
          : 'bg-blue-50 text-blue-700 dark:bg-blue-900/30 dark:text-blue-400'
        }`}>
        {member.role}
      </span>
    </div>
  );

  // --- TAMPILAN SKELETON (JIKA LOADING) ---
  if (isLoading) {
    return (
      <div className="max-w-6xl mx-auto pb-20 px-6 animate-pulse">
        {/* Header Skeleton */}
        <div className="text-center mb-16 pt-10">
          <div className="h-10 bg-gray-200 dark:bg-slate-700 rounded-full w-3/4 md:w-1/2 mx-auto mb-6"></div>
          <div className="h-4 bg-gray-200 dark:bg-slate-700 rounded-full w-full max-w-2xl mx-auto mb-3"></div>
          <div className="h-4 bg-gray-200 dark:bg-slate-700 rounded-full w-2/3 max-w-2xl mx-auto"></div>
        </div>

        {/* Section 1 Grid Skeleton */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mb-24">
           <div className="h-64 bg-gray-200 dark:bg-slate-700 rounded-2xl"></div>
           <div className="space-y-4">
             <div className="h-20 bg-gray-200 dark:bg-slate-700 rounded-xl"></div>
             <div className="h-20 bg-gray-200 dark:bg-slate-700 rounded-xl"></div>
             <div className="h-20 bg-gray-200 dark:bg-slate-700 rounded-xl"></div>
           </div>
        </div>

        {/* Team Skeleton */}
        <div className="flex justify-center gap-6 mb-8">
           <div className="w-64 h-64 bg-gray-200 dark:bg-slate-700 rounded-2xl"></div>
           <div className="w-64 h-64 bg-gray-200 dark:bg-slate-700 rounded-2xl hidden md:block"></div>
        </div>
      </div>
    );
  }

  // --- TAMPILAN ASLI ---
  return (
    <div className="max-w-6xl mx-auto overflow-x-hidden pb-20">
      
      {/* --- HERO HEADER --- */}
      <div className="text-center mb-16 pt-10 px-4">
        <h1 className="text-4xl md:text-5xl font-extrabold text-main mb-6 tracking-tight">
          Tentang PolyChem
        </h1>
        <p className="text-lg text-muted max-w-3xl mx-auto leading-relaxed">
          Aplikasi web untuk memprediksi properti termal polimer — Glass Transition Temperature (Tg) —
          dari struktur molekul SMILES, dengan arsitektur Agentic AI yang menggabungkan validasi kimia dan penalaran model bahasa.
        </p>
      </div>

      {/* --- SECTION 1: LATAR BELAKANG & TUJUAN --- */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 items-start mb-24 px-6">
        
        {/* KARTU KIRI: Latar Belakang */}
        <div className="bg-card p-8 rounded-2xl border border-border shadow-sm hover:shadow-md transition-all h-full">
          <div className="inline-block px-3 py-1 mb-6 text-xs font-semibold tracking-wider text-blue-600 dark:text-blue-400 uppercase bg-blue-50 dark:bg-blue-900/30 rounded-full">
            Latar Belakang
          </div>
          <h2 className="text-2xl font-bold text-main mb-6">Mengapa Riset Ini Penting?</h2>
          <div className="prose dark:prose-invert text-muted leading-relaxed space-y-4 text-justify">
            <p>
              Penemuan material polimer baru secara tradisional (wet lab) memakan waktu lama, biaya tinggi, dan melibatkan proses <em>trial and error</em> yang tidak efisien. Sebagian besar alat prediksi kimia berbasis AI saat ini juga masih berjalan pada antarmuka baris perintah yang sulit diakses peneliti non-programer.
            </p>
            <p>
              PolyChem hadir sebagai aplikasi web berarsitektur <strong>Client-Server</strong> yang mengorkestrasi <strong>Agentic AI</strong>: RDKit memvalidasi struktur SMILES dan mengekstraksi fitur kimia, lalu model bahasa (LLM) melakukan penalaran berbasis fitur tersebut untuk memperkirakan nilai Tg — parameter penting yang menentukan fleksibilitas dan stabilitas termal material polimer.
            </p>
          </div>
        </div>

        {/* KOLOM KANAN: Yang Membedakan PolyChem (Icon Cards) */}
        <div className="grid grid-cols-1 gap-4 h-full">
          {/* Card 1 */}
          <div className="bg-card p-6 rounded-xl border border-border shadow-sm hover:shadow-md transition-shadow flex items-start gap-4">
            <div className="w-12 h-12 rounded-lg bg-blue-100 dark:bg-blue-900/30 flex items-center justify-center text-blue-600 dark:text-blue-400 flex-shrink-0">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" /></svg>
            </div>
            <div>
              <h3 className="font-bold text-main mb-1">Agentic AI Pipeline</h3>
              <p className="text-sm text-muted">RDKit memvalidasi struktur SMILES dan mengekstraksi fitur kimia (Morgan Fingerprint), lalu LLM menalar di atas fitur tersebut untuk meminimalkan halusinasi hasil prediksi.</p>
            </div>
          </div>

          {/* Card 2 */}
          <div className="bg-card p-6 rounded-xl border border-border shadow-sm hover:shadow-md transition-shadow flex items-start gap-4">
            <div className="w-12 h-12 rounded-lg bg-purple-100 dark:bg-purple-900/30 flex items-center justify-center text-purple-600 dark:text-purple-400 flex-shrink-0">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-4.35-4.35M11 19a8 8 0 100-16 8 8 0 000 16z" /></svg>
            </div>
            <div>
              <h3 className="font-bold text-main mb-1">Analisis Kemiripan Molekul</h3>
              <p className="text-sm text-muted">Tanimoto Similarity membandingkan fingerprint senyawa input dengan dataset yang tersedia, menampilkan senyawa serupa sebagai konteks tambahan atas hasil prediksi.</p>
            </div>
          </div>

          {/* Card 3 */}
          <div className="bg-card p-6 rounded-xl border border-border shadow-sm hover:shadow-md transition-shadow flex items-start gap-4">
            <div className="w-12 h-12 rounded-lg bg-green-100 dark:bg-green-900/30 flex items-center justify-center text-green-600 dark:text-green-400 flex-shrink-0">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M20 7l-8-4-8 4m16 0l-8 4m8-4v10l-8 4m0-10L4 7m8 4v10M4 7v10l8 4" /></svg>
            </div>
            <div>
              <h3 className="font-bold text-main mb-1">Arsitektur Client-Server Modular</h3>
              <p className="text-sm text-muted">Frontend (React) dan backend (FastAPI) dipisah secara decoupled, sehingga tiap sisi dapat dikembangkan dan diskalakan secara independen.</p>
            </div>
          </div>
        </div>
      </div>

      {/* --- SECTION 2: TECH STACK (MARQUEE) --- */}
      <div className="mb-24 bg-gray-50 dark:bg-slate-900 py-16 transition-colors duration-300">
        <h2 className="text-2xl font-bold text-main mb-8 text-center">Teknologi di Balik Layar</h2>
        <div className="relative w-full overflow-hidden">
          <div 
            className="flex gap-8 w-max"
            style={{ animation: 'scroll 25s linear infinite' }}
          >
            {[...technologies, ...technologies].map((tech, idx) => (
              <div key={idx} className="w-48 flex-shrink-0 bg-card p-6 rounded-xl border border-border text-center hover:shadow-md transition-all cursor-default flex flex-col items-center justify-center h-40">
                <div className="h-14 w-14 mb-3 flex items-center justify-center">
                  <img src={tech.image} alt={`${tech.name} Logo`} className="h-full w-full object-contain" />
                </div>
                <h3 className="font-bold text-main text-sm">{tech.name}</h3>
                <p className="text-xs text-muted mt-1">{tech.role}</p>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* --- SECTION 3: TIM PENGEMBANG --- */}
      <div className="px-6">
        <div className="text-center mb-12">
          <h2 className="text-3xl font-bold text-main mb-4">Tim Pengembang</h2>
          <p className="text-muted max-w-2xl mx-auto">
            Dikembangkan oleh tim capstone A25-CS110 ASAH By Dicoding 2025
          </p>
        </div>
        
        <div className="flex flex-col items-center gap-8">
          {/* BARIS 1: MACHINE LEARNING (2 Orang) */}
          <div className="flex flex-wrap justify-center gap-6 w-full">
            {mlTeam.map((member, index) => renderMemberCard(member, index))}
          </div>

          {/* BARIS 2: WEB DEVELOPER (*/}
          <div className="flex flex-wrap justify-center gap-6 w-full">
            {webTeam.map((member, index) => renderMemberCard(member, index))}
          </div>
        </div>
      </div>

      {/* --- SECTION 4: PEMBIMBING --- */}
      <div className="px-6 mt-16">
        <div className="text-center mb-12">
          <h2 className="text-3xl font-bold text-main mb-4">Pembimbing</h2>
        </div>

        <div className="flex flex-wrap justify-center gap-6 w-full">
          {advisorTeam.map((member, index) => renderMemberCard(member, index))}
        </div>
      </div>

      {/* STYLE ANIMASI */}
      <style>{`
        @keyframes scroll {
          0% { transform: translateX(0); }
          100% { transform: translateX(-50%); }
        }
      `}</style>

    </div>
  );
}

export default AboutPage;