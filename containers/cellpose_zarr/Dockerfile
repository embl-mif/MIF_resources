# Adapted from the cellposesam-zarr recipe by the EMBL Centre for Bioimage Analysis
# (https://git.embl.org/grp-cba/containers).

# ---- Build arguments ----
ARG PYTHON_VERSION=3.12
ARG CELLPOSE_VERSION=4.1.1
ARG NGFF_ZARR_VERSION=0.34.0
ARG ZARR_VERSION=3.1.6
ARG NUMPY_VERSION=2.4.4

FROM python:${PYTHON_VERSION}-slim

ARG CELLPOSE_VERSION
ARG NGFF_ZARR_VERSION
ARG ZARR_VERSION
ARG NUMPY_VERSION

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PIP_ROOT_USER_ACTION=ignore

# procps provides `ps`, which Nextflow needs to collect task metrics
RUN apt-get update && apt-get install -y \
    build-essential \
    libgl1 \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    numpy==${NUMPY_VERSION} \
    zarr==${ZARR_VERSION} \
    "ngff-zarr[all]==${NGFF_ZARR_VERSION}" \
    aiohttp requests \
    cellpose==${CELLPOSE_VERSION}

CMD ["python"]
