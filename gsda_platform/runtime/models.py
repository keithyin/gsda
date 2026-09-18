"""Runtime model (design doc §20). One long-lived runtime per user in V1."""
from __future__ import annotations

import datetime

from sqlalchemy import DateTime, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column

from gsda_platform.db.database import Base

# §23 — the set of lifecycle states a runtime can be in.
STATUSES = ("CREATING", "RUNNING", "IDLE", "STOPPED", "CRASHED", "FAILED")


def _utcnow() -> datetime.datetime:
    # Naive UTC (SQLite driver drops tzinfo). See auth.models._utcnow.
    return datetime.datetime.utcnow()


class Runtime(Base):
    __tablename__ = "runtimes"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    # UNIQUE enforces the §36 concurrency guard: one runtime per user at the DB level.
    user_id: Mapped[int] = mapped_column(Integer, unique=True, nullable=False, index=True)
    container_id: Mapped[str | None] = mapped_column(Text, nullable=True)
    image: Mapped[str] = mapped_column(Text, nullable=False)
    workspace_path: Mapped[str] = mapped_column(Text, nullable=False)
    results_path: Mapped[str] = mapped_column(Text, nullable=False)
    # dsh state dir (DSH_HOME) — kept out of the doc's schema but needed for §27
    # image/data separation so sessions persist across restarts.
    state_path: Mapped[str] = mapped_column(Text, nullable=False, default="")
    status: Mapped[str] = mapped_column(String(16), nullable=False, default="CREATING")
    crash_restarts: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime(timezone=True), default=_utcnow
    )
    last_active: Mapped[datetime.datetime] = mapped_column(
        DateTime(timezone=True), default=_utcnow
    )
