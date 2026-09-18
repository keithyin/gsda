"""Environment-driven configuration for the Platform Server.

Every knob is an env var (``PLATFORM_*``) with a sane default so the server runs
unmodified in dev, and can be fully configured by the Dockerfile / compose / entrypoint.
"""
from __future__ import annotations

import os
from dataclasses import dataclass


def _env(name: str, default: str | None = None) -> str | None:
    v = os.environ.get(name)
    return v if v not in (None, "") else default


def _int(name: str, default: int) -> int:
    raw = _env(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        raise ValueError(f"{name} must be an integer, got {raw!r}") from None


@dataclass(frozen=True)
class Settings:
    # Storage root for all platform + user data (§25: id-based, not username).
    data_dir: str
    db_path: str

    # The single image, one version (§5, §40).
    image_name: str
    image_tag: str

    # Network.
    platform_port: int
    harness_port: int
    # The public authority users reach (e.g. "analysis.company.internal"). Forwarded
    # verbatim to dsh's --trusted-host and preserved on proxied requests so its
    # browser-trust fence accepts the platform's host (see harness/README.md).
    public_authority: str
    # Bridge network both the platform and user containers join, so the platform can
    # reach a container by IP.
    docker_network: str

    # Runtime lifecycle (§24).
    idle_timeout: int          # seconds; idle > this stops the container (keeps data)
    max_crash_restarts: int    # §37: attempts before a runtime is marked FAILED

    # dsh LLM config passed into the harness container's DSH_HOME/settings.yaml.
    # The key itself is read straight from the process env (ANTHROPIC_API_KEY) so it
    # never lives in Settings or the DB.

    # ------------------------------------------------------------------ paths ----
    @property
    def image(self) -> str:
        return f"{self.image_name}:{self.image_tag}"

    def users_root(self) -> str:
        return os.path.join(self.data_dir, "users")

    def user_dirs(self, user_id: int) -> tuple[str, str, str]:
        """Return (workspace, results, dsh_state) host dirs for a user (§25–27).

        ``dsh_state`` is the per-user ``DSH_HOME`` so dsh sessions/storages survive
        container stop/restart and image upgrades.
        """
        base = os.path.join(self.users_root(), str(user_id))
        return (
            os.path.join(base, "workspace"),
            os.path.join(base, "results"),
            os.path.join(base, "dsh-home"),
        )

    def ensure_user_dirs(self, user_id: int) -> tuple[str, str, str]:
        paths = self.user_dirs(user_id)
        for p in paths:
            os.makedirs(p, exist_ok=True)
        return paths


def load_settings() -> Settings:
    data_dir = _env("PLATFORM_DATA_DIR", "/data")
    return Settings(
        data_dir=data_dir,
        db_path=_env("PLATFORM_DB_PATH", os.path.join(data_dir, "platform.db")),
        image_name=_env("PLATFORM_IMAGE_NAME", "company/analysis-agent"),
        image_tag=_env("PLATFORM_IMAGE_TAG", "1.0.0"),
        platform_port=_int("PLATFORM_PORT", 8000),
        harness_port=_int("PLATFORM_HARNESS_PORT", 3080),
        public_authority=_env("PLATFORM_PUBLIC_AUTHORITY", "analysis.company.internal"),
        docker_network=_env("PLATFORM_DOCKER_NETWORK", "gsda-platform"),
        idle_timeout=_int("PLATFORM_IDLE_TIMEOUT", 1800),
        max_crash_restarts=_int("PLATFORM_MAX_CRASH_RESTARTS", 3),
    )
