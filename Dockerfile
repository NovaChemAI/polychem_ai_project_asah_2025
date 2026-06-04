FROM mambaorg/micromamba:1.5.8

WORKDIR /app

# 1) Install env first for better layer caching
COPY polychem-ai_backend/environment.yml /app/environment.yml
RUN micromamba create -y -n app -f /app/environment.yml && micromamba clean --all --yes

# 2) Copy backend source
COPY polychem-ai_backend/ /app/

# micromamba auto-activate env for RUN/CMD after this
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# 3) Run uvicorn
CMD ["bash", "-lc", "micromamba run -n app uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000} --proxy-headers --forwarded-allow-ips='*'"]
