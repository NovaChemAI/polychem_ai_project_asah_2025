# PolyChem AI
ID Use Case : AC-04
<br>
Use Case    : Novel Chemicals Discovery Agent
<br>
<details>
<summary>📑 Daftar Isi</summary>

- [Latar Belakang](#latar-belakang)
- [Proyek](#proyek)
- [Arsitektur Sistem](#arsitektur-sistem)
- [Model AI yang digunakan](#model-ai-yang-digunakan)
- [Limitasi Sistem](#limitasi-sistem)
- [Fitur](#fitur)
- [Teknologi yang digunakan](#teknologi-yang-digunakan)
- [Instalasi dan Menjalankan Proyek](#instalasi-dan-menjalankan-proyek)
- [Tim Capstone](#tim-capstone)

</details>

## Latar Belakang
Industri petrokimia membutuhkan inovasi senyawa baru, tetapi prosesnya terbilang cukup lambat. Industri ini membutuhkan teknologi yang dapat membantu para ahli kimia dalam efisiensi penemuan senyawa baru. Salah satunya yaitu industri yang memproduksi plastik atau senyawa kimia polymer dalam skala besar. 
Isu-isu terkait kerusakan bumi yang diakibatkan sampah plastik atau senyawa sintetis polimer semakin menguat. Beberapa senyawa polimer pun terindikasi menyebabkan gangguan kesehatan dan berbahaya untuk pemakaian jangka panjang. Oleh karena itu kami membuat proyek ini dengan tujuan agar dapat mempermudah para ahli kimia pada industri produksi plastik dalam menemukan senyawa polimer baru yang dapat diproduksi dengan skala besar, tanpa menyebabkan kerusakan jangka panjang.

## Proyek
<img src="https://drive.google.com/uc?export=view&id=1LKmYSqnHF32he0VMGrZHxw1ArMqiLBeL" width="850">

Dataset : [Data Polimer](https://www.kaggle.com/datasets/linyeping/extra-dataset-with-smilestgpidpolimers-class ) <br>
Link Preprocessing Dataset : [clean dataset](https://colab.research.google.com/drive/1xx0RhQ6OmJE_ZOy46y-f4z2wXTRGIZuc#scrollTo=NinH9YSZiFoI) <br>
Link Model Agentic AI :  [Modelling](https://colab.research.google.com/drive/1IBfzseIefJsAwU7sggxj1tQFYhid5CPj?usp=sharing) <br>

### Arsitektur Sistem
Proyek ini menggunakan arsitektur client–server yang memisahkan antara frontend, backend, dan AI service untuk menjaga skalabilitas dan maintainability.

🔹 Frontend <br>
Frontend dibangun menggunakan:
- React + Vite
- Tailwind CSS
- Firebase Authentication 

Frontend bertugas untuk:
1. Menyediakan antarmuka pengguna (UI)
2. Mengelola interaksi pengguna
3. Mengirim permintaan (request) ke backend
4. Menampilkan hasil respon dari sistem AI

Frontend di-deploy menggunakan Vercel.

🔹 Backend & AI Service
Backend dikembangkan menggunakan FastAPI (Python) dan berfungsi sebagai:
- API server
- Penghubung antara frontend dan sistem AI
- Tempat eksekusi Agentic AI

Backend dan AI service di-deploy menggunakan Koyeb, yang menyediakan lingkungan cloud untuk menjalankan aplikasi Python, pengelolaan environment variable, serta dukungan autoscaling.

```
User → Frontend (React + Vite) → Backend (FastAPI) → Agentic AI → Response
          ^                                                           |
          |-----------------------------------------------------------|
```

### Model AI yang digunakan
Sistem AI pada proyek ini diimplementasikan sebagai custom agentic AI yang dijalankan di backend menggunakan FastAPI. Seluruh alur pengambilan keputusan, pemrosesan data, dan pemanggilan model AI dikendalikan secara manual melalui fungsi Python, tanpa menggunakan framework agent otomatis.

**Komponen utama AI meliputi :** <br><br>
**1️. Large Language Model (LLM) --> Menggunakan Google Generative AI** <br>
Berfungsi untuk:
- Memahami input pengguna
- Memberikan respons berbasis bahasa alami
- Membantu interpretasi hasil analisis molekul

Integrasi LLM dilakukan melalui library google-generativeai. Library LangChain tidak digunakan sebagai agent, melainkan hanya sebagai utility untuk konfigurasi dan manajemen API key.

**2️. Molecular Processing (RDKit)** <br>
Digunakan untuk:
- Parsing struktur molekul (SMILES)
- Pembuatan fingerprint molekul
- Representasi dan visualisasi struktur kimia
- Preprocessing data kimia sebelum analisis
RDKit dijalankan sepenuhnya di backend karena membutuhkan komputasi dan dependensi Python khusus.

**3️. Similarity Measurement (Tanimoto Similarity)**<br>
Berfungsi untuk:
- Mengukur tingkat kemiripan antar molekul
- Memberikan hasil yang lebih relevan untuk senyawa kimia dan polimer
- Perhitungan similarity dilakukan menggunakan fingerprint hasil RDKit.

**4️. Custom Agent Logic (Function-based)**<br>
Alur kerja AI tidak bergantung pada framework agent eksternal, melainkan dibangun secara manual dengan pendekatan berbasis fungsi, seperti:
- Validasi input
- Preprocessing molekul
- Perhitungan similarity
- Pemanggilan LLM
- Postprocessing dan penyusunan respons

Pendekatan ini memberikan kontrol penuh terhadap alur eksekusi dan memudahkan penyesuaian sesuai kebutuhan domain.

**5️. API Key**<br>
- API memakai Google API key, dengan gemini 2.5-flash sebagai model llm
- Google API Key dikelola sepenuhnya di backend

LangChain digunakan hanya sebagai helper untuk integrasi model, bukan sebagai engine agent

### Limitasi Sistem
Seluruh proses kecerdasan buatan pada sistem ini dijalankan di backend, termasuk pemrosesan molekul menggunakan RDKit, perhitungan kemiripan molekul dengan Tanimoto Similarity, serta inferensi menggunakan Large Language Model (LLM). Konsekuensinya, performa sistem sangat bergantung pada kapasitas server dan kestabilan koneksi jaringan, terutama ketika terjadi peningkatan jumlah permintaan.

Sistem ini juga bergantung pada layanan LLM eksternal, yaitu Google Generative AI, sehingga ketersediaan layanan, batas penggunaan API, dan perubahan kebijakan dari penyedia dapat memengaruhi keberlanjutan sistem. Selain itu, pendekatan agentic yang digunakan masih bersifat sederhana karena alur pengambilan keputusan dibangun secara manual tanpa framework agent otomatis, sehingga belum mendukung reasoning multi-step yang adaptif.

Pada aspek analisis kimia, penggunaan Tanimoto Similarity bersifat heuristik dan bergantung pada fingerprint molekul yang dihasilkan oleh RDKit. Hasil analisis tidak mempertimbangkan validasi eksperimental atau sifat fisik molekul secara mendalam, serta sistem terbatas pada format input tertentu seperti SMILES.

### Fitur
- Login page
- Regist page
- Forgot pasword 
- home page 
- statistik 
- library 
- Wiki
- About
- history
- darkmode
- logout
  
### Teknologi yang digunakan
- VSCode
- Python
- Tailwind
- FireBase
- Google Colab
- FastAPI
- Google API
- Vercel
- Koyeb

## Instalasi dan Menjalankan Proyek
Link Deployment : [Klik Disini](https://polychem-ai-project-asah-2025-flame.vercel.app/)

### Front-End
1. Clone repository Github
```
git clone https://github.com/NovaChemAI/polychem_ai_project_asah_2025.git
cd polychem_ai_project_asah_2025/polychem-ai_frontend
```
2. Instalasi depedensi
```
npm install
```
3. Jalankan
```
npm run dev
```
### Back-End
1. Clone repository Github
```
git clone https://github.com/NovaChemAI/polychem_ai_project_asah_2025.git
cd polychem_ai_project_asah_2025/polychem-ai_frontend
```
2. Buat virtual environment
```
python -m venv venv
```
3. Aktifkan
```
Windowa :
venv\Scripts\activate

Mac / Linux :
source venv/bin/activate
```
4. Install dependensi
```
pip install -r requirements.txt
```
5. Jalankan
```
uvicorn main:app --reload
```

## Tim Capstone
ID Team    : A25-CS110 <br>
ID Advisor : A25-ML034 <br>
Advisor    : Yuda Hendriawan Budi Handoko <br>
Anggota tim :

| ID          | Nama                   | Role        |       GitHub      |
|-------------|------------------------|-------------|-------------------|
| M429D5X1913 | Tiara Diansyah Putri   | ML Engineer | @TiaraDian        |
| M156D5X1002 | Loista Amanda Noviar   | ML Engineer | @LoistaAmanda     |
| F284D5Y0483 | Dikky Juliyanto        | Back-End    | @DikkyJuliyanto47 |
| F671D5Y2010 | Zakiul Fata            | Front-End   | @zakiulFata       |

