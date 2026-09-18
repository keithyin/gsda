# Harness contract — DeepSeek Harness (`dsh`)

This directory holds the platform's **launch configuration and contract notes** for the
harness role. It does **not** reimplement the agent — `dsh` (DeepSeek Harness) is the
Agent Runtime (design doc §28). Everything in this file was verified against the installed
`dsh` (v0.1.5-rc.2) rather than assumed.

## Launch

The harness role is started by `docker/entrypoint.sh`:

```bash
cd /workspace
dsh web --host 0.0.0.0 --port "$HARNESS_PORT" --no-open --trusted-host "$PUBLIC_AUTHORITY"
```

with env: `DSH_HOME=/dsh-home`, `DSH_PERMISSION_MODE=workspace-write`.

## Verified facts

| Fact | Value | Source |
|---|---|---|
| Launch | `dsh web` (alias of `dsh --profile web`) | `dsh --help` |
| Flags | `--host`, `--port`, `--no-open`, `--trusted-host <authority...>` | `dsh web --help` |
| Default bind | `127.0.0.1:3080` | `dsh-web-app/cordis.patch.yml` |
| Health endpoint | **none** — no `/health` route exists | grep of the package tree |
| `/api` access | **browser-trust fence** keyed on the request authority | `cordis.patch.yml` |
| Workspace root | `process.cwd()` | `dsh-base/cordis.patch.yml` |
| Write sandbox | `DSH_PERMISSION_MODE` (default `workspace-write`) | `dsh-base/cordis.patch.yml` |
| State home | `DSH_HOME` (default `~/.dsh`): `settings.yaml`, `profiles/`, `sessions/`, `storages/` | filesystem |
| LLM config | `settings.yaml` → `llm-pi-ai.providers`, `agent-default-model`; key via `apiKeyEnv` | `~/.dsh/settings.yaml` |
| Extra env | `DSH_TOOLS_MODE` (native\|ptc\|both) | `dsh-web-app/cordis.patch.yml` |

## Three things that must be handled (and where they are)

1. **No health endpoint.** dsh has no `/health`. The Platform's readiness probe
   (`RuntimeManager` → `DockerBackend.readiness_probe`) is therefore **TCP connect +
   HTTP-200 on `/`**, never a `/health` GET. A container only becomes `RUNNING` once that
   succeeds.

2. **`--trusted-host` is mandatory.** dsh's web app checks the request **Host** against its
   `--trusted-host` value (its browser-trust fence). The browser talks to the *Platform's*
   public host, which reverse-proxies to the container. If the container isn't told that
   same authority, dsh rejects every proxied `/api` call. So:
   - the Platform passes `PLATFORM_PUBLIC_AUTHORITY` into the container as `PUBLIC_AUTHORITY`;
   - the entrypoint feeds it to `dsh web --trusted-host`;
   - the proxy forwards the client's original `Host` **unchanged** (`proxy/routes.py` does
     not rewrite it to the container IP).
   All three must agree. This is the single most likely silent failure in the system.

3. **`--host 0.0.0.0` is required.** The default binds loopback-only, which is unreachable
   from the Platform container. The entrypoint always passes `0.0.0.0`.

## State persistence

`DSH_HOME` is mounted from the user's persistent state dir
(`{DATA_DIR}/users/<id>/dsh-home`) so `profiles/`, `sessions/`, and `storages/` survive
container stop/restart and image upgrades. On first run the Runtime Manager writes a
`settings.yaml` into it (never overwriting one the user has edited) so the harness has an
LLM provider out of the box.

## LLM configuration

The LLM provider is configured in `$DSH_HOME/settings.yaml`, not via dsh flags. The
Runtime Manager seeds it from `ANTHROPIC_BASE_URL` / `PLATFORM_LLM_MODEL` /
`PLATFORM_LLM_PROVIDER`; the key is read from the `ANTHROPIC_API_KEY` process env (passed
through to the container) and is never stored in the platform database.
