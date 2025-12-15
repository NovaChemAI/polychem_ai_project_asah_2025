Dockerfile
FROM mambaorg/micromamba:1.5.8

WORKDIR /app

# Copy env first (cache build)
COPY environment.yml /app/environment.yml
RUN micromamba create -y -n app -f /app/environment.yml && micromamba clean --all --yes

# Copy source
COPY . /app

# Koyeb akan set PORT; pastikan uvicorn bind ke 0.0.0.0
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PYTHONUNBUFFERED=1

CMD ["bash", "-lc", "micromamba run -n app uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}"]