"""Admin routes (design doc §38, §48). Admin-only; ordinary users cannot reach them."""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Request, status
from sqlalchemy import select
from sqlalchemy.orm import Session as SASession

from gsda_platform.auth.models import User
from gsda_platform.deps import get_current_admin, get_db
from gsda_platform.runtime import models as rt_models
from gsda_platform.runtime.manager import RuntimeManager

router = APIRouter(prefix="/api/admin", tags=["admin"])


def _manager(request: Request) -> RuntimeManager:
    return request.app.state.runtime_manager


def _rt_dict(rt: rt_models.Runtime) -> dict:
    return {
        "id": rt.id,
        "user_id": rt.user_id,
        "container_id": rt.container_id,
        "image": rt.image,
        "status": rt.status,
        "crash_restarts": rt.crash_restarts,
        "last_active": rt.last_active.isoformat() if rt.last_active else None,
    }


@router.get("/users")
def list_users(db: SASession = Depends(get_db), _admin: User = Depends(get_current_admin)):
    users = db.scalars(select(User).order_by(User.id)).all()
    return [
        {
            "id": u.id,
            "username": u.username,
            "is_admin": u.is_admin,
            "created_at": u.created_at.isoformat() if u.created_at else None,
        }
        for u in users
    ]


@router.get("/runtimes")
def list_runtimes(db: SASession = Depends(get_db), _admin: User = Depends(get_current_admin)):
    rts = db.scalars(select(rt_models.Runtime)).all()
    return [_rt_dict(rt) for rt in rts]


@router.post("/runtime/{user_id}/restart")
def restart_runtime(user_id: int, request: Request, _admin: User = Depends(get_current_admin)):
    return _rt_dict(_manager(request).restart(user_id).runtime)


@router.post("/runtime/{user_id}/stop")
def stop_runtime(user_id: int, request: Request, _admin: User = Depends(get_current_admin)):
    return _rt_dict(_manager(request).stop(user_id).runtime)


@router.delete("/runtime/{user_id}")
def delete_runtime(user_id: int, request: Request, _admin: User = Depends(get_current_admin)):
    _manager(request).delete(user_id)
    return {"ok": True}
