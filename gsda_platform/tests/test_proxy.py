"""Proxy tests (design doc §32–35, §54).

A real upstream HTTP server is run on a free localhost port; the fake backend's
container IP is pointed at it so the platform's catch-all route forwards to it. This
proves the HTTP proxy path (including streaming) and the authority pass-through that
dsh's trust fence depends on.
"""
from __future__ import annotations

import socket
import threading

import httpx
import pytest
from fastapi import FastAPI, Request


def _free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def _make_upstream() -> tuple[FastAPI, str]:
    upstream = FastAPI()
    seen_host = {}

    @upstream.get("/")
    def root(request: Request):
        seen_host["host"] = request.headers.get("host")
        return {"from": "harness"}

    @upstream.get("/stream")
    def stream():
        import asyncio

        async def gen():
            for i in range(3):
                yield f"chunk{i}\n"
                await asyncio.sleep(0.01)

        from starlette.responses import StreamingResponse

        return StreamingResponse(gen(), media_type="text/plain")

    return upstream, seen_host


@pytest.fixture
def upstream():
    port = _free_port()
    app, seen_host = _make_upstream()
    import uvicorn

    config = uvicorn.Config(app, host="127.0.0.1", port=port, log_level="error")
    server = uvicorn.Server(config)
    t = threading.Thread(target=server.run, daemon=True)
    t.start()
    # wait for readiness
    for _ in range(100):
        try:
            httpx.get(f"http://127.0.0.1:{port}/", timeout=0.2)
            break
        except httpx.HTTPError:
            threading.Event().wait(0.02)
    yield port, seen_host
    server.should_exit = True


def test_proxy_forwards_and_preserves_authority(client, fake_docker, upstream):
    port, seen_host = upstream
    # Log in, create the runtime, then point the fake container at our upstream.
    client.post("/api/auth/register", json={"username": "alice", "password": "password123"})
    client.get("/api/runtime")
    cid = next(iter(fake_docker.containers))
    fake_docker.containers[cid]["ip"] = "127.0.0.1"
    # override the port the proxy targets
    import dataclasses

    client.app.state.settings = dataclasses.replace(client.app.state.settings, harness_port=port)

    r = client.get("/", headers={"host": "analysis.company.internal"})
    assert r.status_code == 200
    assert r.json() == {"from": "harness"}
    # The public Host must survive the proxy (dsh trust fence requirement).
    assert seen_host.get("host") == "analysis.company.internal"


def test_proxy_streams_chunks(client, fake_docker, upstream):
    port, _ = upstream
    client.post("/api/auth/register", json={"username": "alice", "password": "password123"})
    client.get("/api/runtime")
    cid = next(iter(fake_docker.containers))
    fake_docker.containers[cid]["ip"] = "127.0.0.1"
    import dataclasses

    client.app.state.settings = dataclasses.replace(client.app.state.settings, harness_port=port)

    r = client.get("/stream")
    assert r.status_code == 200
    assert r.text == "chunk0\nchunk1\nchunk2\n"


def test_proxy_requires_auth(client):
    r = client.get("/", headers={"host": "analysis.company.internal"})
    assert r.status_code == 401


def test_platform_health_not_shadowed_by_proxy(client):
    """/health must hit the platform's own route, not fall through to the proxy
    (which would require auth). Guards against router-ordering regressions."""
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["ok"] is True
