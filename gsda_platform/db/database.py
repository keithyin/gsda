"""Database setup: SQLAlchemy 2.0 engine/session on SQLite.

V1 uses synchronous SQLAlchemy against a single SQLite file (stdlib sqlite3 driver).
FastAPI route handlers that touch the DB are ``def`` (not ``async def``) so they run in
the threadpool and don't block the event loop; the proxy stays ``async def`` because it
streams. ``aiosqlite`` is left available if we ever go fully async.
"""
from __future__ import annotations

import os

from sqlalchemy import create_engine
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker


class Base(DeclarativeBase):
    pass


def make_engine(db_path: str):
    """Create a SQLite engine. ``check_same_thread=False`` because FastAPI serves the
    sync handlers from a worker threadpool rather than the request thread."""
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    engine = create_engine(
        f"sqlite:///{db_path}",
        connect_args={"check_same_thread": False},
        future=True,
    )
    return engine


def init_db(engine) -> None:
    """Create all tables (idempotent). Import models so their tables are registered."""
    # noqa: F401 — importing for side effects (table registration)
    from gsda_platform.auth import models as _auth_models
    from gsda_platform.runtime import models as _runtime_models  # noqa: F401

    Base.metadata.create_all(engine)


def make_session_factory(engine) -> sessionmaker[Session]:
    return sessionmaker(bind=engine, autoflush=False, expire_on_commit=False)
