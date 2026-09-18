"""Platform Server entrypoint (design doc §7, §47).

The FastAPI app: builds the engine/session-factory at startup, wires the runtime
manager with its Docker backend, and mounts the auth / runtime / admin / proxy routers.
The proxy router is mounted last so it only catches paths that no ``/api/*`` route
handled.

Run locally: ``uvicorn gsda_platform.main:app``.
Inside the image, the platform role runs this on ``PLATFORM_PORT`` (see
``docker/entrypoint.sh``).
"""
from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager

from fastapi import FastAPI

from gsda_platform import __version__
from gsda_platform.config import load_settings
from gsda_platform.db import database

log = logging.getLogger("gsda_platform")


def create_app(settings=None, docker_backend=None, *, start_reaper: bool = False) -> FastAPI:
    """Build the app. ``settings`` and ``docker_backend`` can be injected for tests;
    in production they come from env / the Docker SDK."""
    if settings is None:
        settings = load_settings()

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        engine = database.make_engine(settings.db_path)
        database.init_db(engine)
        app.state.settings = settings
        app.state.engine = engine
        app.state.session_factory = database.make_session_factory(engine)

        if docker_backend is None:
            from gsda_platform.runtime.docker import RealDockerBackend

            backend = RealDockerBackend()
        else:
            backend = docker_backend
        app.state.runtime_manager = _build_manager(app.state.session_factory, settings, backend)

        reaper = None
        if start_reaper:
            reaper = app.state.runtime_manager.start_reaper()
        # The harness reads ANTHROPIC_API_KEY straight from the process env; if it's
        # unset/empty the container boots fine but every LLM call fails at runtime with a
        # confusing auth error. Surface it loudly at startup instead.
        if not os.environ.get("ANTHROPIC_API_KEY"):
            log.warning(
                "ANTHROPIC_API_KEY is not set — harness LLM calls will fail at runtime. "
                "Set it in the environment before serving users."
            )
        log.info("gsda_platform %s up; image=%s", __version__, settings.image)
        try:
            yield
        finally:
            if reaper is not None:
                app.state.runtime_manager.stop_reaper()
            engine.dispose()

    app = FastAPI(title="gsda Platform", version=__version__, lifespan=lifespan)
    # Exposed for observability / tests: whether this app instance will start the
    # idle-reaper loop (§24) at startup.
    app.state.reaper_enabled = start_reaper

    from gsda_platform.auth.routes import router as auth_router
    from gsda_platform.runtime.routes import router as runtime_router
    from gsda_platform.admin.routes import router as admin_router
    from gsda_platform.proxy.routes import router as proxy_router

    app.include_router(auth_router)
    app.include_router(runtime_router)
    app.include_router(admin_router)

    @app.get("/health", tags=["platform"])
    def platform_health():
        # The *platform's* own health endpoint (not the harness's). Registered before
        # the proxy catch-all so it is not swallowed by the ``/{path:path}`` route.
        return {"ok": True, "version": __version__}

    # Last: the proxy's catch-all must not shadow the platform's own routes.
    app.include_router(proxy_router)

    return app


def _build_manager(session_factory, settings, docker_backend):
    from gsda_platform.runtime.manager import RuntimeManager

    return RuntimeManager(session_factory, settings, docker_backend)


# Default app for ``uvicorn gsda_platform.main:app``.
# start_reaper=True so the idle auto-reclaim loop (§24) actually runs in the deployed
# server; the create_app() default stays False so tests/lib consumers opt in explicitly.
app = create_app(start_reaper=True)
