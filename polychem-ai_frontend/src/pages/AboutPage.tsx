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
    { name: "Python AI", role: "Machine Learning", image: pythonLogo },
    { name: "Tailwind CSS", role: "Styling System", image: tailwindLogo },
    { name: "Google Colab", role: "Model Training", image: colabLogo },
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
      <h3 className="font-bold text-lg text-main mb-1">{member.name}</h3>
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
          Tentang PolyChemAI
        </h1>
        <p className="text-lg text-muted max-w-3xl mx-auto leading-relaxed">
          Platform penemuan material cerdas yang menggabungkan kekuatan Artificial Intelligence dengan ilmu polimer untuk masa depan yang lebih berkelanjutan.
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
              Penemuan material polimer baru secara tradisional memakan waktu tahunan dan biaya laboratorium yang sangat besar. Metode <em>trial and error</em> fisik seringkali tidak efisien untuk mengejar kebutuhan material ramah lingkungan yang mendesak.
            </p>
            <p>
              PolyChemAI hadir untuk memotong komputasi tersebut menggunakan <strong>Machine Learning</strong>. Dengan memprediksi properti kimia langsung dari struktur molekul (SMILES), peneliti dapat menyaring ribuan kandidat material hanya dalam hitungan detik.
            </p>
          </div>
        </div>

        {/* KOLOM KANAN: Tujuan Project (Icon Cards) */}
        <div className="grid grid-cols-1 gap-4 h-full">
          {/* Card 1 */}
          <div className="bg-card p-6 rounded-xl border border-border shadow-sm hover:shadow-md transition-shadow flex items-start gap-4">
            <div className="w-12 h-12 rounded-lg bg-blue-100 dark:bg-blue-900/30 flex items-center justify-center text-blue-600 dark:text-blue-400 flex-shrink-0">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" /></svg>
            </div>
            <div>
              <h3 className="font-bold text-main mb-1">Integrasi AI Modern</h3>
              <p className="text-sm text-muted">Menggabungkan model Python canggih dengan antarmuka web yang responsif dan cepat.</p>
            </div>
          </div>

          {/* Card 2 */}
          <div className="bg-card p-6 rounded-xl border border-border shadow-sm hover:shadow-md transition-shadow flex items-start gap-4">
            <div className="w-12 h-12 rounded-lg bg-purple-100 dark:bg-purple-900/30 flex items-center justify-center text-purple-600 dark:text-purple-400 flex-shrink-0">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z" /></svg>
            </div>
            <div>
              <h3 className="font-bold text-main mb-1">User Interface Intuitif</h3>
              <p className="text-sm text-muted">Desain yang ramah pengguna, memudahkan kimiawan tanpa latar belakang coding.</p>
            </div>
          </div>

          {/* Card 3 */}
          <div className="bg-card p-6 rounded-xl border border-border shadow-sm hover:shadow-md transition-shadow flex items-start gap-4">
            <div className="w-12 h-12 rounded-lg bg-green-100 dark:bg-green-900/30 flex items-center justify-center text-green-600 dark:text-green-400 flex-shrink-0">
              <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" /></svg>
            </div>
            <div>
              <h3 className="font-bold text-main mb-1">Kolaborasi Riset</h3>
              <p className="text-sm text-muted">Fitur Library untuk menyimpan dan berbagi riwayat analisis antar peneliti.</p>
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
            Dikembangkan dengan semangat inovasi oleh mahasiswa Universitas Islam Negeri Ar-Raniry sebagai bagian dari Capstone Project.
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