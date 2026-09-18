"""Authentication routes (design doc §48)."""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Response, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session as SASession

from gsda_platform.auth import service
from gsda_platform.auth.models import User
from gsda_platform.deps import clear_session_cookie, cookie, get_current_user, get_db, set_session_cookie

router = APIRouter(prefix="/api/auth", tags=["auth"])


class RegisterIn(BaseModel):
    username: str = Field(min_length=3, max_length=32)
    password: str = Field(min_length=8)


class LoginIn(BaseModel):
    username: str = Field(min_length=1, max_length=64)
    password: str = Field(min_length=1)


@router.post("/register", status_code=status.HTTP_201_CREATED)
def register(body: RegisterIn, response: Response, db: SASession = Depends(get_db)):
    try:
        user = service.register(db, body.username, body.password)
    except service.UsernameTaken:
        raise HTTPException(status.HTTP_409_CONFLICT, "username already taken")
    except ValueError as e:
        raise HTTPException(status.HTTP_422_UNPROCESSABLE_ENTITY, str(e))

    sess = service.create_session(db, user.id)
    set_session_cookie(response, sess.id)
    return {"id": user.id, "username": user.username}


@router.post("/login")
def login(body: LoginIn, response: Response, db: SASession = Depends(get_db)):
    try:
        user = service.authenticate(db, body.username, body.password)
    except service.InvalidCredentials:
        raise HTTPException(status.HTTP_401_UNAUTHORIZED, "invalid credentials")

    sess = service.create_session(db, user.id)
    set_session_cookie(response, sess.id)
    return {"id": user.id, "username": user.username, "is_admin": user.is_admin}


@router.post("/logout")
def logout(
    response: Response,
    cookie_value: str = Depends(cookie),
    db: SASession = Depends(get_db),
):
    # Read the raw cookie directly so a stale/expired session can still be cleared.
    service.delete_session(db, cookie_value)
    clear_session_cookie(response)
    return {"ok": True}


@router.get("/me")
def me(user: User = Depends(get_current_user)):
    return {
        "id": user.id,
        "username": user.username,
        "is_admin": user.is_admin,
        "created_at": user.created_at.isoformat() if user.created_at else None,
    }
