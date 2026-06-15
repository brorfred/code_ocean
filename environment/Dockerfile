# Code Ocean capsule / local podman compatible
#
# Build (from project root):
#   podman build -t abstemp -f environment/Dockerfile .
#
# Run (locally with podman — mount your data dir and an output dir):
#   podman run --rm \
#     -v "$(pwd)/code/src/abstemp/data:/data:ro" \
#     -v "$(pwd)/figs:/results" \
#     abstemp
#
# In Code Ocean: /data is mounted read-only with input data,
#                /results is mounted writable for output.

FROM debian:bookworm-slim

# ── System packages ──────────────────────────────────────────────────────────
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# ── Pixi ─────────────────────────────────────────────────────────────────────
ENV PIXI_HOME=/usr/local/pixi
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_HOME=/usr/local/pixi sh
ENV PATH="/usr/local/pixi/bin:$PATH"

# Store pixi project environments outside /code so they survive Code Ocean's
# /code volume mount at runtime.
RUN mkdir -p /root/.config/pixi \
 && printf 'detached-environments = "/opt/pixi-envs"\n' \
    > /root/.config/pixi/config.toml

# ── Project ───────────────────────────────────────────────────────────────────
COPY code/ /code/

# Code Ocean standard directories
RUN mkdir -p /data /results

# Install the pixi environment into /opt/pixi-envs/ (reproducible via pixi.lock)
RUN cd /code && pixi install --locked

# ── Runtime ───────────────────────────────────────────────────────────────────
# Non-interactive matplotlib backend (no display needed)
ENV MPLBACKEND=Agg
# Disable Python output buffering so progress is visible in container logs
ENV PYTHONUNBUFFERED=1
# Large NetCDF files are provided via /data mount (Code Ocean) or a local bind mount
ENV ABSTEMP_DATA_DIR=/data

WORKDIR /code
ENTRYPOINT ["/bin/bash", "/code/run"]
