#!/usr/bin/env bash
#
# Single-image, two-role entrypoint (design doc §6, §44).
#
#   ROLE=platform  → start the Platform Server (FastAPI/uvicorn)
#   ROLE=harness   → start the user's DeepSeek Harness (dsh web)
#
# One image, two runtimes. The harness role runs in a per-user container whose
# workspace and dsh state are bind-mounted in by the Runtime Manager.
#
set -euo pipefail

ROLE="${ROLE:-}"

# Where the dsh launcher lives (preinstalled into the image; see Dockerfile).
DSH_BIN="${DSH_BIN:-dsh}"

case "$ROLE" in
  platform)
    # The control plane. Binds the public port; proxies into user containers.
    exec python -m uvicorn gsda_platform.main:app \
      --host "${PLATFORM_HOST:-0.0.0.0}" \
      --port "${PLATFORM_PORT:-8000}" \
      --workers 1
    ;;

  harness)
    # The per-user agent runtime. dsh's browser-trust fence validates the request
    # Host, so --trusted-host MUST be the platform's public authority (harness/README.md).
    WORKSPACE="${WORKSPACE_DIR:-/workspace}"
    export DSH_HOME="${DSH_HOME:-/dsh-home}"
    mkdir -p "$WORKSPACE" "$DSH_HOME"
    cd "$WORKSPACE"   # dsh treats process.cwd() as the workspace root

    exec "$DSH_BIN" web \
      --host 0.0.0.0 \
      --port "${HARNESS_PORT:-3080}" \
      --no-open \
      --trusted-host "${PUBLIC_AUTHORITY:?PUBLIC_AUTHORITY is required for the harness role}"
    ;;

  *)
    echo "ERROR: ROLE must be 'platform' or 'harness' (got: '${ROLE}')" >&2
    exit 1
    ;;
esac
