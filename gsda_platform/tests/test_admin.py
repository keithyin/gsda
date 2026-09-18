"""Admin tests (design doc §38, §48)."""
from __future__ import annotations


def _register_admin(client, username="root"):
    client.post("/api/auth/register", json={"username": username, "password": "password123"})
    # Promote the first user to admin directly in the DB (V1 has no admin-bootstrap API).
    from gsda_platform.auth.models import User
    from sqlalchemy import select

    with client.app.state.session_factory() as db:
        u = db.scalar(select(User).where(User.username == username))
        u.is_admin = True
        db.commit()


def test_admin_can_list_users(client):
    client.post("/api/auth/register", json={"username": "alice", "password": "password123"})
    _register_admin(client, "root")
    # Log in as the admin in a fresh client so its cookie is the admin's.
    r = client.post("/api/auth/login", json={"username": "root", "password": "password123"})
    assert r.status_code == 200
    users = client.get("/api/admin/users").json()
    names = {u["username"] for u in users}
    assert "root" in names and "alice" in names


def test_non_admin_cannot_access_admin(client):
    client.post("/api/auth/register", json={"username": "alice", "password": "password123"})
    # alice is not an admin.
    assert client.get("/api/admin/users").status_code == 403


def test_anonymous_cannot_access_admin(client):
    assert client.get("/api/admin/runtimes").status_code == 401


def test_admin_can_see_runtimes(client):
    client.post("/api/auth/register", json={"username": "root", "password": "password123"})
    _register_admin(client, "root")
    # create a runtime as the admin (who is also a normal user)
    client.get("/api/runtime")
    rts = client.get("/api/admin/runtimes").json()
    assert isinstance(rts, list) and len(rts) >= 1
    assert rts[0]["user_id"] >= 1


def test_admin_can_stop_and_delete_runtime(client):
    client.post("/api/auth/register", json={"username": "root", "password": "password123"})
    _register_admin(client, "root")
    rt = client.get("/api/runtime").json()
    uid = rt["user_id"]

    stopped = client.post(f"/api/admin/runtime/{uid}/stop").json()
    assert stopped["status"] == "STOPPED"

    deleted = client.delete(f"/api/admin/runtime/{uid}")
    assert deleted.status_code == 200
    rts = client.get("/api/admin/runtimes").json()
    assert all(r["user_id"] != uid for r in rts)
