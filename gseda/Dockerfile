# Stage 1: builder
FROM ubuntu:24.04 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-dev build-essential \
    libbz2-dev liblz4-dev liblzma-dev libgomp1 pkg-config \
    gcc g++ cmake curl ca-certificates \
    && rm -rf /var/lib/apt/lists/* && \
    curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y --no-install-recommends nodejs && \
    apt-get purge -y curl && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

COPY . /workspace/gseda
WORKDIR /workspace/gseda

# Install Rust with rsproxy mirror for faster downloads
ENV RUSTUP_DIST_SERVER="https://rsproxy.cn" \
    RUSTUP_UPDATE_ROOT="https://rsproxy.cn/rustup"
RUN mkdir -p ~/.cargo && \
    cat > ~/.cargo/config.toml <<'CARGOEOF'
[source.crates-io]
replace-with = 'rsproxy-sparse'

[source.rsproxy]
registry = "https://rsproxy.cn/crates.io-index"

[source.rsproxy-sparse]
registry = "sparse+https://rsproxy.cn/index/"

[registries.rsproxy]
index = "https://rsproxy.cn/crates.io-index"

[net]
git-fetch-with-cli = true
CARGOEOF
RUN curl --proto '=https' --tlsv1.2 -sSf https://rsproxy.cn/rustup-init.sh | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install Claude Code via npm
RUN npm install -g @anthropic-ai/claude-code

# Stage 2: runtime
# FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip libbz2 liblzma libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# COPY --from=builder /workspace/gseda /workspace/gseda
WORKDIR /workspace/gseda

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD ["python3", "-c", "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')"] || exit 1

CMD ["python3", "-m", "uvicorn", "gseda.server.main:app", "--host", "0.0.0.0", "--port", "8000"]
