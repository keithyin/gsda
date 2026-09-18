"""Runtime Manager tests (design doc §19–27, §36–37, §55).

These drive ``RuntimeManager`` directly against a ``FakeDockerBackend`` (no live
Engine). They cover the lifecycle state machine, the §36 concurrency guard, crash
recovery, and the idle reaper.
"""
from __future__ import annotations

import threading

import pytest

from fastapi.testclient import TestClient
from gsda_platform.db import database
from gsda_platform.runtime.docker import FakeDockerBackend
from gsda_platform.runtime.manager import RuntimeManager


@pytest.fixture
def env(settings):
    engine = database.make_engine(settings.db_path)
    database.init_db(engine)
    sf = database.make_session_factory(engine)
    backend = FakeDockerBackend()
    mgr = RuntimeManager(sf, settings, backend)
    return mgr, backend, sf


def test_ensure_runtime_creates_and_runs(env):
    mgr, backend, _ = env
    h = mgr.ensure_runtime(1)
    assert h.status == "RUNNING"
    assert len(backend.containers) == 1


def test_ensure_runtime_is_idempotent(env):
    mgr, backend, _ = env
    first = mgr.ensure_runtime(1)
    second = mgr.ensure_runtime(1)
    assert first.container_id == second.container_id
    assert len(backend.containers) == 1  # no second container


def test_one_runtime_per_user_unique(env):
    mgr, _, sf = env
    mgr.ensure_runtime(1)
    mgr.ensure_runtime(2)
    with sf() as db:
        from gsda_platform.runtime.models import Runtime
        from sqlalchemy import select

        rows = db.scalars(select(Runtime)).all()
    assert len(rows) == 2
    assert {r.user_id for r in rows} == {1, 2}


def test_reuse_running_skips_restart(env):
    mgr, backend, _ = env
    mgr.ensure_runtime(1)
    calls = list(backend.calls)
    mgr.ensure_runtime(1)  # already running → no new start call
    starts_before = calls.count("start:fake0001x")
    assert backend.calls.count("start:fake0001x") == starts_before  # unchanged


def test_stop_then_start_reuses_container(env):
    mgr, backend, _ = env
    h1 = mgr.ensure_runtime(1)
    mgr.stop(1)
    assert mgr.get(1).status == "STOPPED"
    h2 = mgr.start(1)
    assert h2.status == "RUNNING"
    assert h2.container_id == h1.container_id


def test_crash_on_start_marks_crashed(env):
    mgr, backend, _ = env
    backend.fail_starts = 1  # the upcoming start crashes
    h = mgr.ensure_runtime(1)
    assert h.status == "CRASHED"
    assert backend.containers[h.container_id]["state"] == "exited"


def test_readiness_failure_counts_crash_restarts(env):
    mgr, backend, _ = env
    backend.next_readiness = False  # harness never becomes ready
    h = mgr.ensure_runtime(1)
    assert h.status == "CRASHED"
    with mgr._db() as db:
        from gsda_platform.runtime.models import Runtime
        rt = db.get(Runtime, h.runtime.id)
        assert rt.crash_restarts == 1


def test_crash_retries_then_failed(env, settings):
    mgr, backend, _ = env
    backend.always_fail_start = True  # every start crashes the container
    # max_crash_restarts=3 by default: after 4 failed attempts it must be FAILED.
    statuses = [mgr.start(1).status for _ in range(4)]
    assert statuses == ["CRASHED", "CRASHED", "CRASHED", "FAILED"]


def test_missing_container_is_recreated(env):
    mgr, backend, _ = env
    h1 = mgr.ensure_runtime(1)
    # Simulate the container vanishing from the Docker daemon.
    backend.containers.pop(h1.container_id, None)
    h2 = mgr.ensure_runtime(1)
    assert h2.status == "RUNNING"
    assert h2.container_id != h1.container_id
    assert len(backend.containers) == 1  # old one gone, exactly one container exists


def test_idle_reaper_stops_only_idle_runtimes(env, settings):
    import datetime

    mgr, backend, sf = env
    mgr.ensure_runtime(1)
    mgr.ensure_runtime(2)

    # Age runtime 1 beyond the idle timeout; keep runtime 2 fresh.
    with sf() as db:
        from gsda_platform.runtime.models import Runtime
        from sqlalchemy import select

        rt1 = db.scalar(select(Runtime).where(Runtime.user_id == 1))
        rt1.last_active = datetime.datetime.utcnow() - datetime.timedelta(seconds=settings.idle_timeout + 1)
        db.commit()

    stopped = mgr.reap_idle()
    assert stopped == [1]
    assert mgr.get(1).status == "STOPPED"
    assert mgr.get(2).status == "RUNNING"


def test_workspace_paths_are_id_based(settings, env):
    mgr, _, _ = env
    mgr.ensure_runtime(1001)
    h = mgr.get(1001)
    assert "users/1001/workspace" in h.runtime.workspace_path
    assert "users/1001/results" in h.runtime.results_path


def test_concurrent_ensure_creates_one_runtime(env):
    """§36 — N threads racing for the same user must yield exactly one container."""
    mgr, backend, _ = env
    errors = []

    def worker():
        try:
            mgr.ensure_runtime(42)
        except Exception as e:  # pragma: no cover
            errors.append(e)

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors
    assert len(backend.containers) == 1


def test_dsh_home_is_mounted(env):
    mgr, backend, _ = env
    mgr.ensure_runtime(1)
    cid = next(iter(backend.containers))
    mounts = backend.containers[cid]["mounts"]
    assert "/workspace" in mounts.values()
    assert "/results" in mounts.values()
    assert "/dsh-home" in mounts.values()


def test_env_carries_role_and_trust(env):
    mgr, backend, _ = env
    mgr.ensure_runtime(1)
    env_vars = next(iter(backend.containers.values()))["env"]
    assert env_vars["ROLE"] == "harness"
    assert "PUBLIC_AUTHORITY" in env_vars
    assert env_vars["DSH_PERMISSION_MODE"] == "workspace-write"


def test_get_runtime_resumes_stopped(client, authed):
    """§21/§49 — GET /api/runtime is the bring-up trigger: a stopped runtime must come
    back to RUNNING on a plain GET, not be returned stuck in STOPPED."""
    client, _uid = authed
    assert client.get("/api/runtime").json()["status"] == "RUNNING"
    # Idle-stop path: stop it, then a fresh GET must bring it back up.
    assert client.post("/api/runtime/stop").json()["status"] == "STOPPED"
    again = client.get("/api/runtime").json()
    assert again["status"] == "RUNNING"


def test_reaper_enabled_in_production_app(settings, fake_docker):
    """§24 — the deployed server must actually run the idle-reaper loop.

    Two guards:
      * the module-level production `app` (the one `uvicorn gsda_platform.main:app`
        serves) is built with the reaper enabled — this is what was missing before
        (the default was `start_reaper=False`, so containers were never idle-stopped);
      * `create_app(start_reaper=True)` really starts the reaper thread at startup.
    """
    import threading

    import gsda_platform.main as platform_main
    from gsda_platform.main import create_app

    # (a) The production wiring: uvicorn serves `platform_main.app`; its reaper must
    #     be on, otherwise §24 idle auto-reclaim never runs in the deployed server.
    assert platform_main.app.state.reaper_enabled is True

    # (b) The mechanism: enabling it actually spawns the (daemon) reaper thread, and
    #     shutting down stops it. Driven against injected settings/backend, so no real
    #     Docker/DB is touched.
    app = create_app(settings=settings, docker_backend=fake_docker, start_reaper=True)
    with TestClient(app):
        threads = [
            t for t in threading.enumerate()
            if t.name == "runtime-idle-reaper" and t.is_alive()
        ]
        assert len(threads) == 1
