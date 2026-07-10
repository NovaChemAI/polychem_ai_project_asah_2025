import os
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI

load_dotenv()

api_key = os.getenv("GOOGLE_API_KEY") or os.getenv("GEMINI_API_KEY")

print("API KEY TERBACA:", bool(api_key))

if not api_key:
    raise RuntimeError("API key belum terbaca. Cek file .env kamu.")

llm = ChatGoogleGenerativeAI(
    model="gemini-2.5-flash",
    google_api_key=api_key,
    temperature=0.3,
)

response = llm.invoke("Balas satu kata saja: aktif")

print("HASIL GEMINI:", response.content)