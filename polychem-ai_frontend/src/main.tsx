// src/main.tsx
import React from 'react'
import ReactDOM from 'react-dom/client'
import { BrowserRouter } from 'react-router-dom' // <--- 1. PASTIKAN INI DIIMPOR
import App from './App.tsx'
import './index.css'

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    {/* 2. PASTIKAN <App /> DIBUNGKUS <BrowserRouter> SEPERTI INI: */}
    <BrowserRouter>
      <App />
    </BrowserRouter>
  </React.StrictMode>,
)