"""Runtime Manager (design doc §19–27, §36–37).

The core of the platform: maps a logged-in user to exactly one long-lived harness
runtime/container, and drives its lifecycle.

Design points that carry through:

* **One runtime per user** — enforced twice: a per-user ``threading.Lock`` (single
  process, single-thread-of-truth) and the ``runtimes.user_id UNIQUE`` constraint
  (§36). Two concurrent requests for the same user can never spawn two containers.
* **Readiness, not "running"** — a container being ``running`` is not the harness being
  usable (§35). We gate on the TCP+HTTP-200 readiness probe before reporting RUNNING.
* **Image / data separation** (§27) — the image is shared; each user gets private
  host dirs (workspace / results / dsh-home) bind-mounted in, so data survives
  container stop/restart and image upgrades.
"""
from __future__ import annotations

import logging
import os
import threading
import time
from dataclasses import dataclass

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session as SASession

from gsda_platform.config import Settings
from gsda_platform.runtime import docker as dkr
from gsda_platform.runtime import models
from gsda_platform.runtime.models import Runtime

log = logging.getLogger("gsda_platform.runtime")


class RuntimeNotFoundError(Exception):
    """No runtime row exists for the user (e.g. stop/delete before any bring-up)."""


@dataclass
class RuntimeHandle:
    runtime: Runtime
    container_ip: str | None
    container_id: str | None
    status: str


class RuntimeManager:
    def __init__(self, session_factory, settings: Settings, backend: dkr.DockerBackend):
        self._db = session_factory
        self._s = settings
        self._dkr = backend
        self._user_locks: dict[int, threading.Lock] = {}
        self._locks_guard = threading.Lock()
        self._reaper_stop = threading.Event()

    # ------------------------------------------------------------------ locks --
    def _lock_for(self, user_id: int) -> threading.Lock:
        with self._locks_guard:
            if user_id not in self._user_locks:
                self._user_locks[user_id] = threading.Lock()
            return self._user_locks[user_id]

    # ------------------------------------------------------------- container cfg
    def _container_env(self, user_id: int) -> dict[str, str]:
        _, _, state_dir = self._s.ensure_user_dirs(user_id)
        env = {
            "ROLE": "harness",
            "HARNESS_PORT": str(self._s.harness_port),
            "PUBLIC_AUTHORITY": self._s.public_authority,
            "DSH_HOME": "/dsh-home",  # in-container path (mounted from state_dir)
            "DSH_PERMISSION_MODE": "workspace-write",
        }
        # Pass LLM config through; the API key stays in the process env, never the DB.
        for key in ("ANTHROPIC_BASE_URL", "ANTHROPIC_API_KEY", "PLATFORM_LLM_MODEL",
                    "PLATFORM_LLM_PROVIDER"):
            if os.environ.get(key):
                env[key] = os.environ[key]
        return env

    def _container_mounts(self, user_id: int) -> dict[str, str]:
        workspace, results, state_dir = self._s.ensure_user_dirs(user_id)
        return {
            workspace: "/workspace",
            results: "/results",
            state_dir: "/dsh-home",
        }

    def _seed_dsh_home(self, user_id: int) -> None:
        """Write a settings.yaml into the user's DSH_HOME on first run so the harness
        has an LLM provider out of the box. Never overwrites an existing file (V1 keeps
        the user's state across upgrades)."""
        _, _, state_dir = self._s.ensure_user_dirs(user_id)
        settings_file = os.path.join(state_dir, "settings.yaml")
        if os.path.exists(settings_file):
            return
        base_url = os.environ.get("ANTHROPIC_BASE_URL", "http://192.168.3.30:8910")
        model = os.environ.get("PLATFORM_LLM_MODEL", "qwen3.8:27b")
        provider = os.environ.get("PLATFORM_LLM_PROVIDER", "anthropic")
        body = (
            "llm-pi-ai:\n"
            "  providers:\n"
            f"    {provider}:\n"
            f"      baseURL: {base_url}\n"
            "      models:\n"
            f"        - id: {model}\n"
            f"          name: {model}\n"
            "      apiKeyEnv: ANTHROPIC_API_KEY\n"
            "agent-default-model:\n"
            f"  provider: {provider}\n"
            f"  model: {model}\n"
        )
        with open(settings_file, "w") as f:
            f.write(body)

    # ------------------------------------------------------------- lifecycle ----
    def _start_and_wait(self, rt: Runtime) -> RuntimeHandle:
        """Start the container for ``rt`` and wait for harness readiness.

        Increments crash_restarts on failure and returns the resulting status. Never
        raises. Two distinct failure modes:
          * the container process exits on start (crash) — detected via inspect;
          * the container runs but the harness never becomes ready (probe fails).
        Either way it counts toward the §37 retry cap and eventually FAILED.
        """
        with self._db() as db:
            rt = db.get(Runtime, rt.id)
            info = self._dkr.inspect(rt.container_id)
            if info.state == dkr.MISSING:
                # Container vanished; recreate in place (§37 fallback).
                self._recreate_container(db, rt)
                info = self._dkr.inspect(rt.container_id)
            try:
                self._dkr.start(rt.container_id)
            except Exception:
                log.exception("docker start failed for runtime %s", rt.id)
                self._record_crash(rt)
                db.commit()
                return self._handle(rt, None)

            # After start, confirm the container actually stayed up.
            info = self._dkr.inspect(rt.container_id)
            if info.state != dkr.RUNNING:
                # Crashed on start — do NOT probe a dead container.
                self._record_crash(rt)
                db.commit()
                return self._handle(rt, None)

            ok = self._dkr.readiness_probe(info, self._s.harness_port)
            if ok:
                rt.status = "RUNNING"
                rt.last_active = _now()
            else:
                self._record_crash(rt)
            db.commit()
            return self._handle(rt, info.ip)

    def _record_crash(self, rt: Runtime) -> None:
        """Bump the crash counter and transition to CRASHED, or FAILED past the cap."""
        rt.crash_restarts += 1
        rt.status = "FAILED" if rt.crash_restarts > self._s.max_crash_restarts else "CRASHED"

    def _recreate_container(self, db: SASession, rt: Runtime) -> None:
        mounts = self._container_mounts(rt.user_id)
        env = self._container_env(rt.user_id)
        new = self._dkr.create(
            image=rt.image,
            name=f"gsda-rt-{rt.user_id}",
            network=self._s.docker_network,
            env=env,
            mounts=mounts,
        )
        rt.container_id = new.container_id
        db.commit()

    def _fresh(self, db: SASession, runtime_id: int) -> Runtime:
        """Re-read the runtime after a sub-session may have committed changes to it.

        Each ``self._db()`` call is a distinct session, so a commit made inside
        ``_start_and_wait`` is invisible to this session's identity map until we expire
        it. Expiring + re-selecting forces a DB round-trip for the latest state."""
        db.expire_all()
        return db.get(Runtime, runtime_id)

    def ensure_runtime(self, user_id: int) -> RuntimeHandle:
        """Find-or-create the runtime for a user and make sure it's RUNNING (or return
        its current state). Idempotent and concurrency-safe (§21, §22, §36)."""
        with self._lock_for(user_id):
            return self._ensure_runtime_locked(user_id)

    def _ensure_runtime_locked(self, user_id: int) -> RuntimeHandle:
        """The body of ``ensure_runtime``; the caller must already hold the user's lock
        (so internal callers like ``restart`` don't deadlock re-acquiring it)."""
        with self._db() as db:
            rt = self._get_or_create(db, user_id)
            self._materialize(db, rt)
            fresh = self._fresh(db, rt.id)
            ip = self._dkr.inspect(fresh.container_id).ip if fresh.container_id else None
            return self._handle(fresh, ip)

    def _get_or_create(self, db: SASession, user_id: int) -> Runtime:
        rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
        if rt is not None:
            return rt
        workspace, results, state_dir = self._s.ensure_user_dirs(user_id)
        self._seed_dsh_home(user_id)
        rt = Runtime(
            user_id=user_id,
            image=self._s.image,
            workspace_path=workspace,
            results_path=results,
            state_path=state_dir,
            status="CREATING",
        )
        try:
            db.add(rt)
            db.commit()
            db.refresh(rt)
        except IntegrityError:
            # Lost a race with another thread that created the row first.
            db.rollback()
            return db.scalar(select(Runtime).where(Runtime.user_id == user_id))

        mounts = self._container_mounts(user_id)
        env = self._container_env(user_id)
        info = self._dkr.create(
            image=rt.image,
            name=f"gsda-rt-{user_id}",
            network=self._s.docker_network,
            env=env,
            mounts=mounts,
        )
        rt.container_id = info.container_id
        rt.status = "CREATING"
        db.commit()
        return rt

    def _materialize(self, db: SASession, rt: Runtime) -> str:
        """Bring an existing runtime to RUNNING if possible. Returns the resulting
        status string."""
        info = self._dkr.inspect(rt.container_id) if rt.container_id else None
        if info is None or info.state == dkr.MISSING:
            self._recreate_container(db, rt)
            info = self._dkr.inspect(rt.container_id)

        if info.state == dkr.RUNNING:
            rt.status = "RUNNING"
            rt.last_active = _now()
            db.commit()
            return rt.status

        # stopped / exited / created → (re)start and wait for readiness
        handle = self._start_and_wait(rt)
        return handle.status

    # --------------------------------------------------------------- explicit ---
    def start(self, user_id: int) -> RuntimeHandle:
        with self._lock_for(user_id):
            with self._db() as db:
                rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
                if rt is None:
                    rt = self._get_or_create(db, user_id)
                self._materialize(db, rt)
                fresh = self._fresh(db, rt.id)
                ip = self._dkr.inspect(fresh.container_id).ip if fresh.container_id else None
                return self._handle(fresh, ip)

    def stop(self, user_id: int) -> RuntimeHandle:
        with self._lock_for(user_id):
            with self._db() as db:
                rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
                if rt is None:
                    raise RuntimeNotFoundError(user_id)
                self._dkr.stop(rt.container_id)
                rt.status = "STOPPED"
                db.commit()
                return self._handle(db.get(Runtime, rt.id), None)

    def delete(self, user_id: int) -> None:
        with self._lock_for(user_id):
            with self._db() as db:
                rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
                if rt is None:
                    return
                if rt.container_id:
                    try:
                        self._dkr.remove(rt.container_id)
                    except Exception:
                        log.warning("container already gone for runtime %s", rt.id)
                db.delete(rt)
                db.commit()

    def restart(self, user_id: int) -> RuntimeHandle:
        """Force a full recreate+restart of the container (drops the old one)."""
        with self._lock_for(user_id):
            with self._db() as db:
                rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
                if rt is None:
                    # Already holding this user's lock, so delegate to the locked body
                    # (calling ensure_runtime here would deadlock re-acquiring it).
                    return self._ensure_runtime_locked(user_id)
                if rt.container_id:
                    try:
                        self._dkr.remove(rt.container_id)
                    except Exception:
                        pass
                    rt.container_id = None
                rt.crash_restarts = 0
                db.commit()
                self._recreate_container(db, rt)
                self._materialize(db, rt)
                fresh = self._fresh(db, rt.id)
                ip = self._dkr.inspect(fresh.container_id).ip if fresh.container_id else None
                return self._handle(fresh, ip)

    # ---------------------------------------------------------------- helpers ---
    def _handle(self, rt: Runtime, ip: str | None) -> RuntimeHandle:
        return RuntimeHandle(rt, ip, rt.container_id, rt.status)

    def touch(self, user_id: int) -> None:
        """Record activity (called by the proxy) so the idle reaper sees the user is
        active (§24)."""
        with self._db() as db:
            rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
            if rt is not None:
                rt.last_active = _now()
                db.commit()

    def get(self, user_id: int) -> RuntimeHandle | None:
        with self._db() as db:
            rt = db.scalar(select(Runtime).where(Runtime.user_id == user_id))
            if rt is None:
                return None
            info = self._dkr.inspect(rt.container_id) if rt.container_id else None
            ip = info.ip if info else None
            return RuntimeHandle(rt, ip, rt.container_id, rt.status)

    # --------------------------------------------------------------- reaper -----
    def start_reaper(self, poll_interval: float = 60.0) -> threading.Thread:
        """Background idle-stop. Returns the (daemon) thread; the app calls
        ``stop_reaper()`` at shutdown."""
        self._reaper_stop.clear()

        def _loop():
            while not self._reaper_stop.wait(poll_interval):
                try:
                    self.reap_idle()
                except Exception:
                    log.exception("idle reaper pass failed")

        t = threading.Thread(target=_loop, name="runtime-idle-reaper", daemon=True)
        t.start()
        return t

    def stop_reaper(self) -> None:
        self._reaper_stop.set()

    def reap_idle(self) -> list[int]:
        """Stop containers idle longer than the timeout; keep their data (§24).
        Returns the user ids that were stopped."""
        stopped: list[int] = []
        cutoff = _now() - _td(self._s.idle_timeout)
        with self._db() as db:
            rts = db.scalars(
                select(Runtime).where(
                    Runtime.status.in_(["RUNNING", "IDLE"]),
                    Runtime.last_active < cutoff,
                )
            ).all()
            for rt in rts:
                try:
                    self._dkr.stop(rt.container_id)
                except Exception:
                    continue
                rt.status = "STOPPED"
                db.commit()
                stopped.append(rt.user_id)
        return stopped


def _now():
    import datetime

    return datetime.datetime.utcnow()


def _td(seconds: float):
    import datetime

    return datetime.timedelta(seconds=seconds)
