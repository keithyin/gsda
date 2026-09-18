"""Runtime routes (design doc §48)."""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Request, status
from sqlalchemy.orm import Session as SASession

from gsda_platform.deps import get_current_user
from gsda_platform.runtime.manager import RuntimeManager, RuntimeNotFoundError
from gsda_platform.auth.models import User

router = APIRouter(prefix="/api/runtime", tags=["runtime"])


def _manager(request: Request) -> RuntimeManager:
    return request.app.state.runtime_manager


def _to_dict(handle) -> dict:
    return {
        "id": handle.runtime.id,
        "user_id": handle.runtime.user_id,
        "image": handle.runtime.image,
        "status": handle.status,
        "container_id": handle.container_id,
        "workspace_path": handle.runtime.workspace_path,
        "results_path": handle.runtime.results_path,
        "last_active": handle.runtime.last_active.isoformat()
        if handle.runtime.last_active
        else None,
    }


@router.get("")
def get_runtime(request: Request, user: User = Depends(get_current_user)):
    # §21/§49: GET /api/runtime is the bring-up trigger. A STOPPED/CRASHED runtime must
    # be (re)started and readied here, not merely returned in its current state — so we
    # always go through ensure_runtime (idempotent, concurrency-safe) rather than
    # short-circuiting on a merely-existing row.
    handle = _manager(request).ensure_runtime(user.id)
    if handle.status == "FAILED":
        raise HTTPException(status.HTTP_503_SERVICE_UNAVAILABLE, "runtime failed to start")
    return _to_dict(handle)


@router.post("")
def create_runtime(request: Request, user: User = Depends(get_current_user)):
    handle = _manager(request).ensure_runtime(user.id)
    if handle.status == "FAILED":
        raise HTTPException(status.HTTP_503_SERVICE_UNAVAILABLE, "runtime failed to start")
    return _to_dict(handle)


@router.post("/start")
def start_runtime(request: Request, user: User = Depends(get_current_user)):
    handle = _manager(request).start(user.id)
    if handle.status == "FAILED":
        raise HTTPException(status.HTTP_503_SERVICE_UNAVAILABLE, "runtime failed to start")
    return _to_dict(handle)


@router.post("/stop")
def stop_runtime(request: Request, user: User = Depends(get_current_user)):
    try:
        return _to_dict(_manager(request).stop(user.id))
    except RuntimeNotFoundError:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "no runtime to stop")


@router.delete("")
def delete_runtime(request: Request, user: User = Depends(get_current_user)):
    _manager(request).delete(user.id)
    return {"ok": True}
