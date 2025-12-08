import reactLogo from '../assets/TechLogo/LogoReact.png';
import firebaseLogo from '../assets/TechLogo/LogoFirebase.png';
import pythonLogo from '../assets/TechLogo/LogoPhyton.png'; 
import tailwindLogo from '../assets/TechLogo/LogoTailwind.png';
import colabLogo from '../assets/TechLogo/LogoColab.png';
import vscodeLogo from '../assets/TechLogo/LogoVsCode.png';

// ▼ IMPORT GAMBAR ILUSTRASI TIM ▼
// Perhatikan: Nama file menggunakan spasi, jadi pastikan path-nya tepat sesuai screenshot
import userMale from '../assets/DevImages/ilustration user.jpg';
import userFemale from '../assets/DevImages/ilustration userw.png';

function AboutPage() {
  // Data Teknologi
  const technologies = [
    { name: "React + Vite", role: "Frontend Modern", image: reactLogo },
    { name: "Firebase", role: "Auth & Database", image: firebaseLogo },
    { name: "Python AI", role: "Machine Learning", image: pythonLogo },
    { name: "Tailwind CSS", role: "Styling System", image: tailwindLogo },
    { name: "Google Colab", role: "Model Training", image: colabLogo },
    { name: "VS Code", role: "Code Editor", image: vscodeLogo },
  ];

  // Data Tim (Updated dengan Image)
  const teamMembers = [
    {
      name: "Tiara Diansyah Putri",
      role: "Machine Learning",
      image: userFemale, // Pakai userw.png
      color: "green"
    },
    {
      name: "Loista Amanda Noviar",
      role: "Machine Learning",
      image: userFemale, // Pakai userw.png
      color: "green"
    },
    {
      name: "Zakiul Fata",
      role: "Front-End & Back-End",
      image: userMale, // Pakai user.jpg
      color: "blue"
    },
    {
      name: "M Ihsan Sadik",
      role: "Front-End & Back-End",
      image: userMale, // Pakai user.jpg
      color: "blue"
    },
    {
      name: "Dikky Juliyanto",
      role: "Front-End & Back-End",
      image: userMale, // Pakai user.jpg
      color: "blue"
    },
  ];

  return (
    <div className="max-w-5xl mx-auto overflow-x-hidden">
      
      {/* --- HERO HEADER --- */}
      <div className="text-center mb-16 pt-8 px-4">
        <h1 className="text-4xl font-extrabold text-gray-900 mb-4">Tentang NovaChemAI</h1>
        <p className="text-xl text-gray-600 max-w-3xl mx-auto leading-relaxed">
          Platform penemuan material berbasis kecerdasan buatan yang dirancang untuk mengakselerasi riset polimer dan plastik ramah lingkungan.
        </p>
      </div>

      {/* --- SECTION 1: LATAR BELAKANG --- */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-12 items-start mb-20 px-4">
        <div>
          <h2 className="text-2xl font-bold text-gray-900 mb-4 flex items-center">
            <span className="bg-blue-100 text-blue-600 w-8 h-8 rounded-lg flex items-center justify-center mr-3 text-sm">01</span>
            Latar Belakang Masalah
          </h2>
          <div className="prose text-gray-600 leading-relaxed space-y-4">
            <p>
              Penemuan material polimer baru secara tradisional memakan waktu tahunan dan biaya yang sangat besar di laboratorium.
            </p>
            <p>
              Kami melihat peluang untuk memotong komputasi tersebut menggunakan <strong>Machine Learning</strong>. Dengan memprediksi properti kimia dari struktur molekul (SMILES), peneliti dapat menyaring ribuan kandidat hanya dalam hitungan detik.
            </p>
          </div>
        </div>

        <div className="bg-white p-8 rounded-2xl border border-gray-200 shadow-sm">
          <h3 className="font-bold text-gray-900 mb-6 border-b border-gray-100 pb-2">Tujuan Capstone Project</h3>
          <ul className="space-y-4">
            <li className="flex items-start">
              <div className="bg-green-100 text-green-700 p-1 rounded-full mr-3 mt-0.5">✓</div>
              <span className="text-gray-700">Mengintegrasikan Model AI Python dengan Web App Modern.</span>
            </li>
            <li className="flex items-start">
              <div className="bg-green-100 text-green-700 p-1 rounded-full mr-3 mt-0.5">✓</div>
              <span className="text-gray-700">Menyediakan antarmuka (UI) yang mudah digunakan.</span>
            </li>
            <li className="flex items-start">
              <div className="bg-green-100 text-green-700 p-1 rounded-full mr-3 mt-0.5">✓</div>
              <span className="text-gray-700">Menyimpan riwayat analisis untuk kolaborasi tim riset.</span>
            </li>
          </ul>
        </div>
      </div>

      {/* --- SECTION 2: TECH STACK (MARQUEE) --- */}
      <div className="mb-24">
        <h2 className="text-2xl font-bold text-gray-900 mb-8 text-center">Teknologi Yang Digunakan</h2>
        
        {/* Container Marquee */}
        <div className="relative w-full overflow-hidden py-4">
          <div 
            className="flex gap-8 w-max"
            style={{
              animation: 'scroll 20s linear infinite', 
            }}
          >
            {/* Render Item (Double) */}
            {[...technologies, ...technologies].map((tech, idx) => (
              <div 
                key={idx} 
                className="w-48 flex-shrink-0 bg-white p-6 rounded-xl border border-gray-200 text-center hover:shadow-md transition-all cursor-default flex flex-col items-center justify-center"
              >
                <div className="h-16 w-16 mb-4 flex items-center justify-center">
                  <img src={tech.image} alt={`${tech.name} Logo`} className="h-full w-full object-contain" />
                </div>
                <h3 className="font-bold text-gray-800 text-sm">{tech.name}</h3>
                <p className="text-xs text-gray-500 mt-1">{tech.role}</p>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* --- SECTION 3: TIM PENGEMBANG (FOTO ILUSTRASI) --- */}
      <div className="bg-white rounded-3xl p-12 border border-gray-200 shadow-sm text-center mb-20 mx-4">
        <h2 className="text-3xl font-bold text-gray-900 mb-4">Tim Pengembang</h2>
        <p className="text-gray-600 mb-12 max-w-2xl mx-auto">
          Proyek ini dikembangkan oleh mahasiswa Universitas Islam Negeri Ar-Raniry sebagai bagian dari Capstone Project.
        </p>
        
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-8 justify-center">
          {teamMembers.map((member, index) => (
            <div key={index} className="bg-gray-50 p-6 rounded-2xl border border-gray-200 hover:bg-white hover:shadow-md transition-all group">
              
              {/* Foto Profil Bulat */}
              <div className={`w-24 h-24 rounded-full mx-auto mb-4 p-1 shadow-sm ring-4 ring-white overflow-hidden bg-white
                ${member.color === 'green' ? 'group-hover:ring-green-100' : 'group-hover:ring-blue-100'} transition-all`}>
                <img 
                  src={member.image} 
                  alt={member.name} 
                  className="w-full h-full object-cover rounded-full"
                />
              </div>
              
              <h3 className="font-bold text-lg text-gray-900">{member.name}</h3>
              <p className={`text-sm font-medium mb-2 ${member.color === 'green' ? 'text-green-700' : 'text-blue-700'}`}>
                {member.role}
              </p>
            </div>
          ))}
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