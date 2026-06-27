# ============================================================
# Pipeline de Variantes Germinativas — Dockerfile
# ============================================================
# Multi-stage build para imagem otimizada.
#
# Build:
#   docker build -t pipeline-variantes:2.0.0 -f docker/Dockerfile .
#
# Run:
#   docker run -v ./samples:/data/samples \
#              -v ./references:/data/references \
#              -v ./results:/data/results \
#              pipeline-variantes:2.0.0
# ============================================================

# ---- Stage 1: Build ----
FROM condaforge/mambaforge:latest AS builder

LABEL maintainer="Pipeline de Variantes Germinativas"
LABEL version="2.0.0"
LABEL description="Pipeline de chamada de variantes germinativas GATK4 Best Practices"

# Copia definição do ambiente
COPY environment.yaml /tmp/environment.yaml

# Cria ambiente Conda otimizado
RUN mamba env create -f /tmp/environment.yaml -n pipeline \
    && mamba clean --all --yes

# ---- Stage 2: Runtime ----
FROM debian:bookworm-slim AS runtime

LABEL maintainer="Pipeline de Variantes Germinativas"
LABEL version="2.0.0"

# Dependências mínimas do runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        sqlite3 \
        procps \
        tini \
    && rm -rf /var/lib/apt/lists/*

# Copia o ambiente Conda do builder
COPY --from=builder /opt/conda/envs/pipeline /opt/conda/envs/pipeline

# Adiciona o ambiente ao PATH
ENV PATH="/opt/conda/envs/pipeline/bin:$PATH"
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Cria usuário não-root
RUN groupadd -r bioinfo && useradd -r -g bioinfo -m -s /bin/bash bioinfo

# Diretórios de trabalho
RUN mkdir -p /data/samples /data/references /data/results \
             /data/config /data/logs /data/database \
    && chown -R bioinfo:bioinfo /data

# Copia o código do pipeline
COPY --chown=bioinfo:bioinfo . /opt/pipeline

WORKDIR /opt/pipeline

# Troca para usuário não-root
USER bioinfo

# Variáveis de ambiente para override
ENV PIPELINE_CONFIG="/data/config/config.yaml"
ENV PIPELINE_SAMPLES_DIR="/data/samples"
ENV PIPELINE_RESULTS_DIR="/data/results"
ENV PIPELINE_REFERENCES_DIR="/data/references"

# Entrypoint com tini para signal handling correto
ENTRYPOINT ["tini", "--"]

# Comando padrão: pipeline completo
CMD ["snakemake", "--cores", "all", \
     "--configfile", "/data/config/config.yaml", \
     "--keep-going", "--printshellcmds", "--rerun-incomplete"]

# ---- Healthcheck ----
HEALTHCHECK --interval=60s --timeout=10s --start-period=5s --retries=3 \
    CMD snakemake --version || exit 1
