import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react-swc'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    headers: {
      // Header ini mengizinkan aplikasi Anda berkomunikasi dengan Popup Google
      "Cross-Origin-Opener-Policy": "same-origin-allow-popups",
    },
  },
})