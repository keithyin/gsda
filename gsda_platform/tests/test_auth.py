"""Auth tests (design doc §12–18, §55)."""
from __future__ import annotations


def test_register_and_me(client):
    r = client.post("/api/auth/register", json={"username": "Alice", "password": "password123"})
    assert r.status_code == 201
    me = client.get("/api/auth/me").json()
    assert me["username"] == "alice"  # stored case-folded (§15)
    assert me["is_admin"] is False


def test_username_case_fold_is_duplicate(client):
    client.post("/api/auth/register", json={"username": "Alice", "password": "password123"})
    dup = client.post("/api/auth/register", json={"username": "ALICE", "password": "password123"})
    assert dup.status_code == 409


def test_duplicate_username_rejected(client):
    client.post("/api/auth/register", json={"username": "bob", "password": "password123"})
    assert client.post("/api/auth/register", json={"username": "bob", "password": "anotherpass"}).status_code == 409


def test_short_username_rejected(client):
    r = client.post("/api/auth/register", json={"username": "ab", "password": "password123"})
    assert r.status_code == 422


def test_short_password_rejected(client):
    r = client.post("/api/auth/register", json={"username": "carol", "password": "short"})
    assert r.status_code == 422


def test_password_not_stored_in_plaintext(client):
    """§16 — the DB must hold a hash, never the raw password."""
    client.post("/api/auth/register", json={"username": "dave", "password": "password123"})
    from gsda_platform.auth.models import User
    from sqlalchemy import select

    factory = client.app.state.session_factory
    with factory() as db:
        user = db.scalar(select(User).where(User.username == "dave"))
        assert user.password_hash != "password123"
        assert user.password_hash.startswith("$argon2")


def test_login_success_and_failure(client):
    client.post("/api/auth/register", json={"username": "erin", "password": "password123"})
    assert client.post("/api/auth/login", json={"username": "erin", "password": "password123"}).status_code == 200
    assert client.post("/api/auth/login", json={"username": "erin", "password": "wrongpass"}).status_code == 401
    assert client.post("/api/auth/login", json={"username": "nope", "password": "password123"}).status_code == 401


def test_session_requires_cookie(client):
    assert client.get("/api/auth/me").status_code == 401


def test_logout_invalidates_session(client):
    client.post("/api/auth/register", json={"username": "frank", "password": "password123"})
    assert client.get("/api/auth/me").status_code == 200
    assert client.post("/api/auth/logout").status_code == 200
    assert client.get("/api/auth/me").status_code == 401


def test_cookie_is_httponly(client):
    r = client.post("/api/auth/register", json={"username": "gina", "password": "password123"})
    set_cookie = r.headers.get("set-cookie", "")
    assert "gsda_session=" in set_cookie
    assert "httponly" in set_cookie.lower()
    assert "samesite=lax" in set_cookie.lower()
