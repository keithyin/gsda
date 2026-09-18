"""Harness reverse proxy (design doc §32–35).

Every request the browser makes that is NOT one of the platform's own ``/api/*``
endpoints is forwarded to the authenticated user's harness container, on the dsh web
port. Three transports must all work (§34):

* **HTTP** (including the dsh gateway's ``/api`` — a different namespace from the
  platform's ``/api/auth`` etc.)
* **WebSocket** (dsh's real-time transport)
* **Streaming / SSE** (chunked responses)

The critical subtlety: dsh validates the request **authority** (Host header) against its
``--trusted-host`` value (its browser-trust fence). We therefore forward the request's
original public Host *unchanged* — we do NOT rewrite it to the container's IP:port, or
dsh would reject every proxied ``/api`` call. See harness/README.md.
"""
from __future__ import annotations

import asyncio

import httpx
from fastapi import APIRouter, Depends, HTTPException, Request, WebSocket, WebSocketDisconnect, status
from fastapi.responses import StreamingResponse

from gsda_platform.auth.models import User
from gsda_platform.deps import get_current_user
from gsda_platform.runtime.manager import RuntimeHandle, RuntimeManager

router = APIRouter(tags=["proxy"])

# Host headers we should NOT rewrite: keep whatever the client sent (the public
# authority). We drop hop-by-hop headers on the way in.
_HOP_HEADERS = {
    "connection", "keep-alive", "proxy-authenticate", "proxy-authorization",
    "te", "trailer", "transfer-encoding", "upgrade", "proxy-connection",
}


def _manager(request) -> RuntimeManager:
    return request.app.state.runtime_manager


def _target_ip(handle) -> str:
    if not handle.container_ip:
        raise HTTPException(status.HTTP_503_SERVICE_UNAVAILABLE, "runtime not running")
    return handle.container_ip


async def _ensure_running(request: Request, user: User) -> RuntimeHandle:
    """Make sure the user's harness is RUNNING and return its handle.

    The manager calls (DB + Docker daemon, and ``ensure_runtime`` which can block on a
    TCP+HTTP readiness probe for seconds to minutes) are all synchronous, so they run
    in a worker thread — otherwise a cold-start request would stall the ASGI event
    loop and freeze every other user's concurrent traffic for that duration.
    """
    mgr = _manager(request)

    def _resolve():
        handle = mgr.get(user.id)
        if handle is None or handle.status not in ("RUNNING", "IDLE"):
            handle = mgr.ensure_runtime(user.id)
        if handle.status != "RUNNING":
            raise HTTPException(status.HTTP_503_SERVICE_UNAVAILABLE, f"runtime {handle.status}")
        mgr.touch(user.id)
        return handle

    return await asyncio.to_thread(_resolve)


def _filter_headers(headers) -> dict:
    out = {}
    for k, v in headers.items():
        if k.lower() in _HOP_HEADERS:
            continue
        out[k] = v
    return out


@router.api_route(
    "/{path:path}",
    methods=["GET", "POST", "PUT", "DELETE", "PATCH", "OPTIONS", "HEAD"],
    include_in_schema=False,
)
async def http_proxy(request: Request, user: User = Depends(get_current_user)):
    handle = await _ensure_running(request, user)
    host = _target_ip(handle)
    port = request.app.state.settings.harness_port
    url = f"http://{host}:{port}/{request.url.path.lstrip('/')}"

    headers = _filter_headers(request.headers)
    # Keep the public Host so dsh's trust fence accepts it.
    headers["host"] = request.headers.get("host", request.app.state.settings.public_authority)

    body = await request.body()

    # No read/write/pool timeout (SSE/streaming must not cut off and the proxy may
    # wait on a busy upstream), but a bounded connect timeout so an unreachable/filtered
    # harness IP can't hang the request and hold a client open indefinitely.
    client = httpx.AsyncClient(
        timeout=httpx.Timeout(connect=10.0, read=None, write=None, pool=None)
    )
    req = client.build_request(
        request.method,
        url,
        headers=headers,
        content=body,
        params=request.query_params,
    )
    resp = await client.send(req, stream=True)

    async def stream():
        try:
            async for chunk in resp.aiter_bytes():
                yield chunk
        finally:
            await resp.aclose()
            await client.aclose()

    # Pass through dsh's status + headers, dropping hop-by-hop and any content-length
    # that httpx may have already consumed while streaming.
    resp_headers = {
        k: v for k, v in resp.headers.items() if k.lower() not in _HOP_HEADERS
        and k.lower() != "content-length"
    }
    return StreamingResponse(
        stream(),
        status_code=resp.status_code,
        headers=resp_headers,
        media_type=resp.headers.get("content-type"),
    )


@router.websocket("/{path:path}")
async def ws_proxy(websocket: WebSocket, user: User = Depends(get_current_user)):
    """Full-duplex WebSocket bridge to the user's dsh endpoint.

    Relays frames both directions concurrently. If either end closes, the other is
    closed too. On close, propagates the code/reason.
    """
    handle = await _ensure_running(websocket, user)
    host = _target_ip(handle)
    port = websocket.app.state.settings.harness_port

    await websocket.accept()

    query = websocket.url.query
    path = websocket.url.path.lstrip("/")
    url = f"ws://{host}:{port}/{path}" + (f"?{query}" if query else "")

    # Forward any client-supplied subprotocols dsh might have negotiated.
    subprotocols = websocket.scope.get("subprotocols") or None

    from websockets.asyncio.client import connect as ws_connect

    try:
        async with ws_connect(url, subprotocols=subprotocols) as upstream:
            async def client_to_upstream():
                try:
                    while True:
                        msg = await websocket.receive()
                        if msg.get("type") == "websocket.disconnect":
                            await upstream.close(code=msg.get("code", 1000))
                            return
                        if "text" in msg and msg["text"] is not None:
                            await upstream.send(msg["text"])
                        elif "bytes" in msg and msg["bytes"] is not None:
                            await upstream.send(msg["bytes"])
                except WebSocketDisconnect:
                    try:
                        await upstream.close()
                    except Exception:
                        pass

            async def upstream_to_client():
                try:
                    async for frame in upstream:
                        if isinstance(frame, (bytes, bytearray)):
                            await websocket.send_bytes(bytes(frame))
                        else:
                            await websocket.send_text(frame)
                    # upstream sent a close frame and finished iterating
                    await websocket.close()
                except Exception:
                    pass

            await asyncio.gather(client_to_upstream(), upstream_to_client())
    except WebSocketDisconnect:
        pass
    except Exception:
        # upstream unreachable / refused — let the client see a clean close
        try:
            await websocket.close(code=1011)
        except Exception:
            pass
