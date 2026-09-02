FROM mambaorg/micromamba:1.5 as builder

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml

RUN micromamba create -y -p /opt/conda_env -f /tmp/env.yml && \
    micromamba clean --all --yes

FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1

COPY --from=builder /opt/conda_env /opt/conda_env

ENV PATH="/opt/conda_env/bin:/app/resources/usr/bin:${PATH}"
ENV MODULE_LIB="/app/resources/lib"

WORKDIR /app

COPY resources/ /app/resources

RUN chmod +x /app/resources/usr/bin/extract_features.py

# Run as non-root
RUN useradd -m pipeline && \
    chown -R pipeline /app
USER pipeline

CMD ["extract_features.py", "--help"]
