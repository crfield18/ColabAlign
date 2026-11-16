# syntax=docker/dockerfile:1
#
# ColabAlign container built from the exact linux-64 reference lock
# (locks/conda-linux-64.lock: Python 3.14.10, USalign, BeEM, TreeCluster).
#
# Build from the repository root:
#   docker build -f containers/Dockerfile -t colabalign:1.2.0 .
#
# Run on the current directory (inputs and outputs live in /data):
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/data" colabalign:1.2.0 \
#       -i structures -o results -c 8 -t 0.2 0.3
#
# The lock is linux-64 only, so the image is always linux/amd64. On Apple
# Silicon it runs under emulation.

ARG MICROMAMBA_IMAGE=mambaorg/micromamba:2.10.0-debian13-slim@sha256:e1476182cd40eb1f7c0cfae4adec9bfd78f3d7f05c61368306abfdc8939c1198
# hadolint ignore=DL3029
FROM --platform=linux/amd64 ${MICROMAMBA_IMAGE}

LABEL org.opencontainers.image.title="ColabAlign" \
      org.opencontainers.image.description="All-vs-all US-align protein structure comparison, structural dendrograms and clustering." \
      org.opencontainers.image.version="1.2.0" \
      org.opencontainers.image.source="https://github.com/crfield18/ColabAlign" \
      org.opencontainers.image.licenses="MIT"

# Exact, hash-checked environment from the reference lock.
COPY --chown=$MAMBA_USER:$MAMBA_USER locks/conda-linux-64.lock /tmp/conda-linux-64.lock
RUN micromamba install --yes --name base --file /tmp/conda-linux-64.lock \
    && micromamba clean --all --yes \
    && rm /tmp/conda-linux-64.lock

# Install ColabAlign itself without touching the locked dependencies.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md LICENSE MANIFEST.in /tmp/colabalign-src/
COPY --chown=$MAMBA_USER:$MAMBA_USER colabalign /tmp/colabalign-src/colabalign
RUN python -m pip install --no-deps --no-build-isolation --no-cache-dir /tmp/colabalign-src \
    && rm -rf /tmp/colabalign-src \
    && python -m pip check \
    && colabalign --version \
    && USalign -h > /dev/null \
    && command -v BeEM > /dev/null \
    && TreeCluster.py -h > /dev/null

# Ignore any user site-packages from a bind-mounted home directory.
ENV PYTHONNOUSERSITE=1

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "colabalign"]
CMD ["--help"]
