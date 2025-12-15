FROM mambaorg/micromamba:1.5.8

WORKDIR /app

# 1) Install env dulu biar caching build bagus
COPY environment.yml /app/environment.yml
RUN micromamba create -y -n app -f /app/environment.yml && micromamba clean --all --yes

# 2) Copy source (pastikan file permission aman)
COPY . /app

# micromamba auto-activate env untuk RUN/CMD setelah ini
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PYTHONUNBUFFERED=1

# Optional tapi bagus untuk PaaS
ENV PYTHONDONTWRITEBYTECODE=1

# 3) Run uvicorn
# - proxy-headers: cocok di behind reverse proxy (Koyeb)
# - forwarded-allow-ips='*': biar gak nge-block forwarded headers
CMD ["bash", "-lc", "micromamba run -n app uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000} --proxy-headers --forwarded-allow-ips='*'"]