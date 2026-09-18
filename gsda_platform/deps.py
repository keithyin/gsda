"""Shared FastAPI dependencies: current DB session and authenticated user.

The app object (built in ``main.py``) stores its engine/session-factory in
``app.state``; these deps pull from ``request.app`` so tests can mount a throwaway
engine without touching globals.
"""
from __future__ import annotations

from typing import Iterator

from fastapi import Depends, HTTPException, Request, status
from fastapi.security import APIKeyCookie
from sqlalchemy.orm import Session as SASession

from gsda_platform.auth import service as auth_service
from gsda_platform.auth.models import User

COOKIE_NAME = "gsda_session"

cookie = APIKeyCookie(name=COOKIE_NAME, auto_error=False)


def get_db(request: Request) -> Iterator[SASession]:
    factory = request.app.state.session_factory
    with factory() as db:
        yield db


def get_current_user(
    cookie_value: str = Depends(cookie),
    db: SASession = Depends(get_db),
) -> User:
    user = auth_service.get_user_by_session(db, cookie_value)
    if user is None:
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, "not authenticated")
    return user


def get_current_admin(user: User = Depends(get_current_user)) -> User:
    if not user.is_admin:
        raise HTTPException(status.HTTP_403_FORBIDDEN, "admin required")
    return user


def set_session_cookie(response, session_id: str) -> None:
    response.set_cookie(
        COOKIE_NAME,
        session_id,
        httponly=True,
        samesite="lax",
        secure=False,  # V1 is internal HTTP; flip to True behind TLS later (§17)
    )


def clear_session_cookie(response) -> None:
    response.delete_cookie(COOKIE_NAME)
