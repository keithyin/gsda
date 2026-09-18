# gsda Platform — single image, two roles (design doc §5–6, §27, §29–31).
#
# One image (`company/analysis-agent:<tag>`) carries both the Platform Server and the
# DeepSeek Harness. `docker/entrypoint.sh` picks which to run from `ROLE`.
#
# Build:
#   docker build -t company/analysis-agent:1.0.0 .
#
# NOTE: not buildable/tested on this machine (no Docker Engine). Review before first
# real build. Node + dsh are preinstalled so containers start without network access.

ARG PYTHON_VERSION=3.10
FROM python:${PYTHON_VERSION}-slim

# Node for dsh (DeepSeek Harness). Pin a known-good LTS.
ARG NODE_VERSION=22
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates build-essential \
    && curl -fsSL https://deb.nodesource.com/setup_${NODE_VERSION}.x | bash - \
    && apt-get install -y nodejs \
    && apt-get remove -y build-essential \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/analysis

# --- Python deps -------------------------------------------------------------
# The bioinformatics stack (gseda/gsetl tooling) plus the platform's web stack.
COPY requirements.txt ./
# `docker` + `argon2-cffi` are the two platform-only additions (see requirements.txt).
RUN pip install --no-cache-dir -r requirements.txt

# --- Analysis capabilities (design §29–31) ----------------------------------
# Mirror the repo's analysis tree under /opt/analysis, PRESERVING the original
# sub-paths (.claude/skills, scripts, third_party). Several skill files hardcode
# absolute host paths (e.g. /root/projects/gsda/.claude/skills/.../run.py and
# /root/projects/gsda/third_party/gseda/...); because the layout is preserved, one
# substitution rewrites all of them to valid in-image paths.
COPY .claude/   /opt/analysis/.claude/
COPY scripts/   /opt/analysis/scripts/
COPY third_party/ /opt/analysis/third_party/

# Rewrite the hardcoded host paths so every skill resolves inside the image.
RUN set -eux; \
    # §29 names the packaged skills dir /opt/analysis/skills — expose it via a
    # symlink so dsh can be pointed there while the real tree keeps .claude/skills.
    ln -s /opt/analysis/.claude/skills /opt/analysis/skills; \
    grep -rl '/root/projects/gsda' /opt/analysis/.claude/skills 2>/dev/null | \
    xargs -r sed -i 's#/root/projects/gsda#/opt/analysis#g' || true; \
    # Fail the build if any host path survived (would silently break every skill).
    ! grep -rq '/root/projects/gsda' /opt/analysis/.claude/skills

# --- Platform server + harness contract -------------------------------------
COPY gsda_platform/   ./gsda_platform/
COPY docker/entrypoint.sh /opt/analysis/docker/entrypoint.sh
RUN chmod +x /opt/analysis/docker/entrypoint.sh

# --- Preinstall dsh so containers don't need npx at boot --------------------
# `npm i -g` caches into the image; the entrypoint calls `dsh` directly.
# Pinned (was `@latest`, which made rebuilds non-deterministic). Bump deliberately.
RUN npm install -g @deepseek-ai/dsh@0.1.5-rc.2

# --- User data dirs (bind-mounted at runtime, §25–27) -----------------------
# These exist so a misconfigured mount still has a target; the Runtime Manager
# mounts real per-user dirs over them.
RUN mkdir -p /workspace /results /dsh-home

ENV PATH="/opt/analysis:${PATH}"
ENV DSH_HOME=/dsh-home
EXPOSE 8000 3080

ENTRYPOINT ["/opt/analysis/docker/entrypoint.sh"]
