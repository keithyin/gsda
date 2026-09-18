"""Shared fixtures for the platform test suite."""
from __future__ import annotations

import dataclasses

import pytest

from gsda_platform.config import load_settings
from gsda_platform.main import create_app
from gsda_platform.runtime.docker import FakeDockerBackend
from fastapi.testclient import TestClient


@pytest.fixture
def settings(tmp_path):
    s = load_settings()
    return dataclasses.replace(
        s,
        data_dir=str(tmp_path / "data"),
        db_path=str(tmp_path / "platform.db"),
        idle_timeout=1800,
        max_crash_restarts=3,
    )


@pytest.fixture
def fake_docker():
    return FakeDockerBackend()


@pytest.fixture
def app(settings, fake_docker):
    """A fully-wired app against a throwaway SQLite DB and a fake Docker backend."""
    return create_app(settings=settings, docker_backend=fake_docker, start_reaper=False)


@pytest.fixture
def client(app):
    with TestClient(app) as c:
        yield c


@pytest.fixture
def authed(client):
    """A client already logged in as 'alice' (returns (client, user_id))."""
    client.post("/api/auth/register", json={"username": "alice", "password": "password123"})
    return client, 1
