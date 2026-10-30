# Adapted from the multicut_stitching recipe by the EMBL Centre for Bioimage Analysis
# (https://git.embl.org/grp-cba/containers).

# ---- Build arguments ----
ARG PYTHON_VERSION=3.12
ARG NUMPY_VERSION=2.4.4
ARG ZARR_VERSION=3.1.6
ARG NIFTY_VERSION=1.2.4
ARG BIOIMAGE_CPP_VERSION=0.8.0

FROM mambaorg/micromamba:1.5.10

ARG PYTHON_VERSION
ARG NUMPY_VERSION
ARG ZARR_VERSION
ARG NIFTY_VERSION
ARG BIOIMAGE_CPP_VERSION

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# procps provides `ps`, which Nextflow needs to collect task metrics.
# The base image runs as $MAMBA_USER, so switch to root for apt.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# nifty and bioimage-cpp are conda-forge-only (the "nifty" PyPI package is an
# unrelated probabilistic-programming library), so install everything via micromamba.
RUN micromamba install -y -n base -c conda-forge \
        python=${PYTHON_VERSION} \
        numpy=${NUMPY_VERSION} \
        zarr=${ZARR_VERSION} \
        nifty=${NIFTY_VERSION} \
        bioimage-cpp=${BIOIMAGE_CPP_VERSION} && \
    micromamba clean -a -y && \
    find /opt/conda -follow -type f -name '*.a'   -delete || true && \
    find /opt/conda -follow -type f -name '*.pyc' -delete || true

CMD ["python"]
