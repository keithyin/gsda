"""Docker backend abstraction (the mock seam).

``RuntimeManager`` talks to containers only through this protocol, so its full
lifecycle logic (create / reuse / crash-restart / idle-stop / readiness) is unit-testable
without a live Docker Engine — tests supply a ``FakeDockerBackend``. ``RealDockerBackend``
wraps the ``docker`` SDK for production.
"""
from __future__ import annotations

import socket
import time
import urllib.request
from dataclasses import dataclass, field
from typing import Protocol

# Container states as Docker reports them.
RUNNING = "running"
STOPPED = "stopped"
EXITED = "exited"
MISSING = "missing"


@dataclass
class ContainerInfo:
    container_id: str
    state: str  # running / stopped / exited / missing
    short_id: str = ""
    ip: str | None = None


def _probe_until_ready(host: str, port: int, *, timeout: float = 60.0) -> bool:
    """Block until ``host:port`` answers TCP + HTTP-200 on ``/``, else time out.

    dsh has NO /health endpoint — readiness is a TCP connect plus a 200 on ``/``
    (see harness/README.md). Shared by the Protocol default and the real backend so
    the loop lives in one place.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with socket.create_connection((host, port), timeout=2) as s:
                s.sendall(b"GET / HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n\r\n")
                first = s.recv(128)
            if first.startswith(b"HTTP/1.1 200"):
                return True
        except OSError:
            pass
        time.sleep(0.5)
    return False


class DockerBackend(Protocol):
    def create(
        self,
        *,
        image: str,
        name: str,
        network: str,
        env: dict[str, str],
        mounts: dict[str, str],  # host_path -> container_path
        hostname: str | None = None,
    ) -> ContainerInfo:
        """Create (not start) a container. Returns its info; id may be a short id."""

    def start(self, container_id: str) -> None: ...

    def stop(self, container_id: str) -> None: ...

    def remove(self, container_id: str) -> None: ...

    def inspect(self, container_id: str) -> ContainerInfo: ...

    def readiness_probe(self, info: ContainerInfo, port: int, *, timeout: float = 60.0) -> bool:
        """Block until the harness answers TCP+HTTP-200 on ``port``, else timeout.

        dsh has NO /health endpoint — readiness is a TCP connect plus a 200 on ``/``
        (see harness/README.md). The default impl here is real; tests override it.
        """
        return _probe_until_ready(info.ip or "127.0.0.1", port, timeout=timeout)


# ------------------------------------------------------------------ real -------
class RealDockerBackend:
    """Production backend over the ``docker`` SDK. Requires a reachable Engine (V1:
    ``/var/run/docker.sock`` mounted into the platform container, §10).

    The client is created **lazily** on first use: the Platform Server must be able to
    boot (auth, admin, DB) even before the Docker Engine is reachable; the connection
    is only needed when a runtime is actually created/started.
    """

    def __init__(self, docker_host: str | None = None):
        self._docker_host = docker_host
        self._client = None

    @property
    def _docker(self):
        if self._client is None:
            import docker

            self._client = (
                docker.from_env()
                if self._docker_host is None
                else docker.DockerClient(base_url=self._docker_host)
            )
        return self._client

    def create(self, *, image, name, network, env, mounts, hostname=None) -> ContainerInfo:
        host_cfg = {c: {"bind": h, "mode": "rw"} for h, c in mounts.items()}
        container = self._docker.containers.create(
            image,
            command=None,  # ENTRYPOINT/ROLE come from the image; env carries ROLE
            name=name,
            hostname=hostname,
            environment=env,
            network=network,
            host_config=self._docker.api.create_host_config(binds=host_cfg),
            detach=True,
        )
        return self._info(container)

    def start(self, container_id: str) -> None:
        self._get(container_id).start()

    def stop(self, container_id: str) -> None:
        self._get(container_id).stop(timeout=10)

    def remove(self, container_id: str) -> None:
        c = self._get(container_id)
        try:
            c.stop(timeout=10)
        except Exception:
            pass
        c.remove(force=True)

    def inspect(self, container_id: str) -> ContainerInfo:
        return self._info(self._get(container_id))

    def readiness_probe(self, info: ContainerInfo, port: int, *, timeout: float = 60.0) -> bool:
        return _probe_until_ready(info.ip or "127.0.0.1", port, timeout=timeout)

    # -- helpers --------------------------------------------------------------
    def _get(self, container_id: str):
        # Accept full or short id; inspect is lenient.
        try:
            return self._docker.containers.get(container_id)
        except Exception:
            raise KeyError(container_id)

    def _info(self, container) -> ContainerInfo:
        state = container.attrs["State"]["State"]
        ip = None
        n = container.attrs.get("NetworkSettings", {}).get("Networks", {})
        for net in n.values():
            if net.get("IPAddress"):
                ip = net["IPAddress"]
                break
        cid = container.id
        return ContainerInfo(
            container_id=cid,
            state=state,
            short_id=cid[:12],
            ip=ip,
        )


# ------------------------------------------------------------------ fake -------
class FakeDockerBackend:
    """In-memory backend for tests. Records calls, models state transitions, and lets
    tests script readiness outcomes and crash behaviour."""

    def __init__(self):
        self.containers: dict[str, dict] = {}
        self._counter = 0
        self.calls: list[str] = []
        self.next_readiness = True  # what readiness_probe returns
        self.fail_starts = 0        # number of upcoming start() calls that crash
        self.always_fail_start = False  # once set, every start crashes the container

    def create(self, *, image, name, network, env, mounts, hostname=None) -> ContainerInfo:
        self._counter += 1
        cid = f"fake{self._counter:04d}x"
        self.calls.append(f"create:{name}")
        self.containers[cid] = {
            "image": image,
            "name": name,
            "network": network,
            "env": env,
            "mounts": mounts,
            "state": "created",
            "ip": f"172.18.0.{self._counter}",
        }
        return ContainerInfo(cid, "created", cid[:12], f"172.18.0.{self._counter}")

    def start(self, container_id: str) -> None:
        self.calls.append(f"start:{container_id}")
        c = self._containers[container_id]
        if self.fail_starts > 0:
            self.fail_starts -= 1
            c["state"] = "exited"  # crash on start
            return
        if self.always_fail_start:
            c["state"] = "exited"
            return
        c["state"] = "running"

    def stop(self, container_id: str) -> None:
        self.calls.append(f"stop:{container_id}")
        self._containers[container_id]["state"] = "stopped"

    def remove(self, container_id: str) -> None:
        self.calls.append(f"remove:{container_id}")
        self.containers.pop(container_id, None)

    def inspect(self, container_id: str) -> ContainerInfo:
        c = self._containers.get(container_id)
        if c is None:
            return ContainerInfo(container_id, "missing")
        state = c["state"] if c["state"] != "created" else "stopped"
        return ContainerInfo(container_id, state, container_id[:12], c.get("ip"))

    def readiness_probe(self, info: ContainerInfo, port: int, *, timeout: float = 60.0) -> bool:
        self.calls.append("readiness")
        return self.next_readiness

    @property
    def _containers(self) -> dict[str, dict]:
        return self.containers
