"""Authentication service (design doc §12–18).

Pure logic, no HTTP. Routes (§48) call these. Raises small domain exceptions the route
layer translates into status codes:

* ``UsernameTaken``      → 409
* ``InvalidCredentials`` → 401
"""
from __future__ import annotations

import datetime
import re
import secrets
import uuid

from argon2 import PasswordHasher
from argon2.exceptions import VerifyMismatchError
from sqlalchemy import select
from sqlalchemy.orm import Session as SASession

from gsda_platform.auth.models import Session, User

SESSION_TTL = datetime.timedelta(days=7)

_USERNAME_RE = re.compile(r"^[a-zA-Z0-9_-]{3,32}$")
_MIN_PASSWORD = 8

_hasher = PasswordHasher()  # argon2id by default


class UsernameTaken(Exception):
    pass


class InvalidCredentials(Exception):
    pass


def _utcnow() -> datetime.datetime:
    return datetime.datetime.utcnow()


def _fold(username: str) -> str:
    """§15 — case-fold so Alice / alice / ALICE are one user."""
    return username.strip().lower()


def validate_username(username: str) -> str:
    if not _USERNAME_RE.match(username or ""):
        raise ValueError("username must be 3-32 chars of [a-zA-Z0-9_-]")
    return _fold(username)


def _validate_password(password: str) -> None:
    if len(password or "") < _MIN_PASSWORD:
        raise ValueError("password must be at least 8 characters")


def register(db: SASession, username: str, password: str, *, is_admin: bool = False) -> User:
    uname = validate_username(username)
    _validate_password(password)

    existing = db.scalar(select(User).where(User.username == uname))
    if existing is not None:
        raise UsernameTaken(uname)

    user = User(username=uname, password_hash=_hasher.hash(password), is_admin=is_admin)
    db.add(user)
    db.commit()
    db.refresh(user)
    return user


def verify_password(stored_hash: str, password: str) -> bool:
    try:
        return _hasher.verify(stored_hash, password)
    except VerifyMismatchError:
        return False


def create_session(db: SASession, user_id: int, *, ttl: datetime.timedelta = SESSION_TTL) -> Session:
    now = _utcnow()
    sess = Session(
        id=secrets.token_urlsafe(32),
        user_id=user_id,
        created_at=now,
        expires_at=now + ttl,
        last_active=now,
    )
    db.add(sess)
    db.commit()
    return sess


def get_valid_session(db: SASession, session_id: str | None) -> Session | None:
    """Return the session if it exists and is not expired, else None."""
    if not session_id:
        return None
    sess = db.get(Session, session_id)
    if sess is None:
        return None
    if sess.expires_at < _utcnow():
        db.delete(sess)
        db.commit()
        return None
    return sess


def get_user_by_session(db: SASession, session_id: str | None) -> User | None:
    sess = get_valid_session(db, session_id)
    if sess is None:
        return None
    return db.get(User, sess.user_id)


def authenticate(db: SASession, username: str, password: str) -> User:
    """Login: verify credentials. Raises InvalidCredentials on any failure."""
    try:
        uname = validate_username(username)
    except ValueError:
        raise InvalidCredentials()
    user = db.scalar(select(User).where(User.username == uname))
    if user is None or not verify_password(user.password_hash, password):
        raise InvalidCredentials()
    return user


def delete_session(db: SASession, session_id: str | None) -> None:
    if not session_id:
        return
    sess = db.get(Session, session_id)
    if sess is not None:
        db.delete(sess)
        db.commit()
