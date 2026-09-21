# --- Stage 1: The Builder ---
FROM condaforge/miniforge3 AS builder

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /build

# 1. Install system build tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# 2. Create the conda environment
COPY environment.yml .
RUN mamba env create -p /opt/conda_env -f environment.yml && \
    mamba clean -afy

# 3. Install PyTorch and Pymif (which handles almost all other dependencies)
# We install torch first to ensure the +cu121 version is prioritized
RUN /opt/conda_env/bin/pip install --no-cache-dir \
    torch torchvision --extra-index-url https://download.pytorch.org/whl/cu121 \
    cellpose \
    git+https://github.com/grinic/pymif.git@v2026.6.0

# --- Stage 2: The Runner (Final Image) ---
FROM debian:bookworm-slim

ENV PATH="/opt/conda_env/bin:$PATH"
ENV CONDA_PREFIX="/opt/conda_env"
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# Install runtime libraries (essential for OpenCV and image processing in pymif)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libfontconfig1 \
    libglx0 \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libgomp1 \
    procps \
    && rm -rf /var/lib/apt/lists/*

    # Copy the environment from the builder
COPY --from=builder /opt/conda_env /opt/conda_env

# Create non-root user
RUN useradd -m nextflow_user
USER nextflow_user

CMD ["python3"]
