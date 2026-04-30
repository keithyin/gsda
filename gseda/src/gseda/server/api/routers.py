"""API Routers - Define all API endpoints for CLI tools"""

import asyncio
import json
import os
import sys
import tempfile
import uuid
from typing import List, Optional, Dict
from fastapi import APIRouter, HTTPException, status, UploadFile, File
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from gseda.server.config import settings
from gseda.server.core.runners import CLIRunner, ARGUMENT_SCHEMAS, get_cli_executor, FileRegistry
from fastapi.responses import FileResponse as FastAPIFileResponse
from gseda.server.core.schema import (
    ToolExecutionRequest,
    ToolExecutionResponse,
    ToolListResponse,
    ToolInfoResponse,
    APIResponse,
)
from gseda.server.core.files import FileManager
from gseda.server.api.error_analysis_using_cc import run_claude_code


def print_log(message: str):
    """Print log message to console with timestamp"""
    import datetime
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


# Create router
api_router = APIRouter(prefix="/tools", tags=["tools"])


# ============================================================================
# File Management Endpoints
# ============================================================================


class RemoteFetchRequest(BaseModel):
    """Request for fetching remote file via SCP"""

    remote_path: str
    ssh_password: str
    ssh_user: Optional[str] = "root"


class RemoteFetchResponse(BaseModel):
    """Response for remote file fetch"""

    local_path: str
    filename: str


class UploadResponse(BaseModel):
    """Response for file upload"""

    local_path: str
    filename: str
    size: int


@api_router.post(
    "/files/remote-fetch",
    response_model=APIResponse,
    summary="Fetch file from remote server via SCP",
    description="Download a file from a remote server using SSH/SCP",
)
async def remote_fetch(request: RemoteFetchRequest) -> APIResponse:
    """Fetch a file from a remote server via SCP."""
    try:
        file_manager = FileManager()
        local_path = file_manager.fetch_remote_file(
            remote_path=request.remote_path,
            ssh_password=request.ssh_password,
            ssh_user=request.ssh_user,
        )

        return APIResponse(
            success=True,
            data=RemoteFetchResponse(
                local_path=local_path,
                filename=os.path.basename(local_path),
            ).dict(),
        )
    except ValueError as e:
        return APIResponse(
            success=False,
            error="Remote fetch failed",
            message=str(e),
        )
    except Exception as e:
        return APIResponse(
            success=False,
            error="Remote fetch failed",
            message=f"Unexpected error: {str(e)}",
        )


@api_router.post(
    "/files/upload",
    response_model=APIResponse,
    summary="Upload a file from client",
    description="Upload a file from the client to the server",
)
async def upload_file(file: UploadFile = File(...)) -> APIResponse:
    """Upload a file from client."""
    try:
        file_manager = FileManager()

        # Generate temp file path
        ext = os.path.splitext(file.filename)[1]
        unique_name = f"gseda_{uuid.uuid4().hex}{ext}"
        local_path = os.path.join(tempfile.gettempdir(), unique_name)

        # Save uploaded file
        with open(local_path, "wb") as f:
            content = await file.read()
            f.write(content)

        # Get file size
        file_size = os.path.getsize(local_path)

        # Track for cleanup
        FileManager._uploaded_files.append(local_path)

        return APIResponse(
            success=True,
            data=UploadResponse(
                local_path=local_path,
                filename=file.filename,
                size=file_size,
            ).dict(),
        )
    except Exception as e:
        return APIResponse(
            success=False,
            error="Upload failed",
            message=str(e),
        )


# ============================================================================
# Tool File Download Endpoints
# ============================================================================


class ToolFileDownloadRequest(BaseModel):
    """Request for downloading tool output file"""
    filename: str


@api_router.get(
    "/files/download/{filename}",
    summary="Download tool output file",
    description="Download a file output from a previously executed CLI tool",
)
async def download_tool_file(filename: str) -> FastAPIFileResponse:
    """Download a tool output file by filename."""
    import os as _os

    # Path traversal protection
    if ".." in filename or filename.startswith("/"):
        raise HTTPException(status_code=400, detail="Invalid filename")

    from gseda.server.core.runners import _load_file_index

    index = _load_file_index()
    entry = None
    for e in index:
        entry_fname = e.get("filename") or e.get("name", "")
        if entry_fname == filename and _os.path.exists(e.get("path", "")):
            entry = e
            break

    if not entry:
        raise HTTPException(
            status_code=404, detail=f"File '{filename}' not found")

    return FastAPIFileResponse(
        path=entry["path"],
        filename=filename,
        media_type="application/octet-stream",
    )


# ============================================================================
# Tool Execution Endpoints
# ============================================================================


@api_router.get(
    "/list",
    response_model=APIResponse,
    summary="List all available CLI tools",
    description="Returns a list of all available CLI tools that can be executed via the API",
)
async def list_tools() -> APIResponse:
    """Get list of all available tools"""
    tools = [
        {"name": name, "module_path": module_path}
        for name, module_path in settings.TOOLS_CONFIG
    ]
    # Add binary tools to the list
    for name, bin_path in settings.TOOL_BINARIES.items():
        tools.append({"name": name, "module_path": bin_path})
    return APIResponse(success=True, data={"tools": tools})


@api_router.get(
    "/{tool_name}/info",
    response_model=APIResponse,
    summary="Get tool information",
    description="Returns metadata and argument schema for a specific tool",
)
async def get_tool_info(tool_name: str) -> APIResponse:
    """Get information about a specific CLI tool"""
    # Find the module path for the tool
    module_path = None
    for name, module in settings.TOOLS_CONFIG:
        if name == tool_name:
            module_path = module
            break
    # Check binary tools
    if not module_path and tool_name in settings.TOOL_BINARIES:
        module_path = settings.TOOL_BINARIES[tool_name]

    if not module_path:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Tool '{tool_name}' not found",
        )

    # Get argument schema if available
    args_schema = ARGUMENT_SCHEMAS.get(tool_name, {"arguments": []})

    tool_info = ToolInfoResponse(
        name=tool_name,
        module_path=module_path,
        description=f"CLI tool: {tool_name}",
        arguments=args_schema.get("arguments", []),
    )

    return APIResponse(success=True, data=tool_info.dict())


# ============================================================================
# AI Analysis Endpoint
# ============================================================================

import asyncio
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor

# In-memory job store for async AI analysis
# key: job_id, value: { "status": "running"|"done"|"error", "result": str, "created_at": float }
AI_ANALYSIS_JOBS: Dict[str, dict] = {}
AI_ANALYSIS_TTL = 300  # jobs expire after 5 minutes
_analysis_thread_pool = ThreadPoolExecutor(max_workers=4)


def _cleanup_expired_jobs():
    """Remove jobs older than TTL."""
    now = time.time()
    expired = [k for k, v in AI_ANALYSIS_JOBS.items() if now - v["created_at"] > AI_ANALYSIS_TTL]
    for k in expired:
        del AI_ANALYSIS_JOBS[k]


def _run_analysis_in_thread(job_id: str, tool_name: str, stderr: str, stdout: str):
    """Background thread that runs the AI analysis for a job."""
    try:
        analysis = f"[工具: {tool_name}]\n\nSTDERR:\n{stderr}"
        if stdout:
            analysis += f"\n\nSTDOUT:\n{stdout}"
        analysis_res = run_claude_code(analysis)
        AI_ANALYSIS_JOBS[job_id] = {
            "status": "done",
            "result": analysis_res or "(AI 未返回结果)",
            "created_at": AI_ANALYSIS_JOBS[job_id]["created_at"],
        }
    except Exception as e:
        AI_ANALYSIS_JOBS[job_id] = {
            "status": "error",
            "result": f"AI 分析失败: {str(e)}",
            "created_at": AI_ANALYSIS_JOBS[job_id]["created_at"],
        }


def _start_analysis_async(job_id: str, tool_name: str, stderr: str, stdout: str):
    """Schedule analysis in the thread pool (non-blocking)."""
    _analysis_thread_pool.submit(_run_analysis_in_thread, job_id, tool_name, stderr, stdout)


class AIAnalyzeRequest(BaseModel):
    """Request for AI analysis of CLI tool stderr"""
    stderr: str
    stdout: Optional[str] = ""
    tool_name: Optional[str] = ""


class AIJobResponse(BaseModel):
    """Response containing a job ID"""
    job_id: str


class AIAnalysisPollResponse(BaseModel):
    """Response for polling AI analysis status"""
    status: str  # "running", "done", "error"
    result: Optional[str] = None


@api_router.post(
    "/{tool_name}/ai-analyze",
    response_model=APIResponse,
    summary="Start AI analysis of CLI error",
    description="Create an AI analysis job for the given stderr/stdout and return a job_id for polling.",
)
async def ai_analyze_start(
    tool_name: str,
    request: AIAnalyzeRequest,
) -> APIResponse:
    """Start async AI analysis of CLI tool stderr. Returns job_id for polling."""
    try:
        # Normalize stderr input: prepend tool name context
        stderr_text = request.stderr or ""
        if request.tool_name:
            stderr_text = f"[工具: {request.tool_name}]\n{stderr_text}"
        else:
            stderr_text = f"[工具: {tool_name}]\n{stderr_text}"

        if request.stdout:
            stderr_text += f"\n\nSTDOUT:\n{request.stdout}"

        # Clean up expired jobs periodically
        _cleanup_expired_jobs()

        # Create job
        job_id = uuid.uuid4().hex[:12]
        AI_ANALYSIS_JOBS[job_id] = {
            "status": "running",
            "result": None,
            "created_at": time.time(),
        }

        # Start analysis in background thread (non-blocking)
        _start_analysis_async(job_id, tool_name, request.stderr or "", request.stdout or "")

        return APIResponse(
            success=True,
            data=AIJobResponse(job_id=job_id).dict(),
        )
    except Exception as e:
        return APIResponse(
            success=False,
            error="Failed to start AI analysis",
            message=str(e),
        )


@api_router.get(
    "/{tool_name}/ai-analyze/{job_id}",
    response_model=APIResponse,
    summary="Poll AI analysis status and result",
    description="Poll the status of an AI analysis job. Returns the result when analysis completes.",
)
async def ai_analyze_poll(
    tool_name: str,
    job_id: str,
) -> APIResponse:
    """Poll the status and result of an AI analysis job."""
    try:
        # Clean up expired jobs on each poll
        _cleanup_expired_jobs()

        job = AI_ANALYSIS_JOBS.get(job_id)
        if not job:
            return APIResponse(
                success=False,
                error="Job not found",
                message=f"Job '{job_id}' not found or has expired.",
            )

        return APIResponse(
            success=True,
            data=AIAnalysisPollResponse(
                status=job["status"],
                result=job["result"],
            ).dict(),
        )
    except Exception as e:
        return APIResponse(
            success=False,
            error="Polling failed",
            message=str(e),
        )


# ============================================================================
# Tool File Download Endpoints
# ============================================================================


def _is_remote_path(path: str) -> bool:
    """
    Check if a path is a remote SCP path.

    SCP format: user@host:/path or host:/path

    Args:
        path: File path to check

    Returns:
        True if path appears to be a remote SCP path
    """
    # Empty or local paths
    if not path or path.startswith('/') or path.startswith('\\') or path.startswith('~'):
        return False

    # Check for user@host format
    if "@" in path and ":" in path:
        return True

    # Check for host:/path format (contains : with nothing before it or non-drive content)
    # Windows paths: C:\ or C:/ (letter followed by colon)
    # Remote paths: host:/path (word followed by colon, no drive letter)
    if ":" in path:
        colon_idx = path.index(":")
        before_colon = path[:colon_idx]
        # If nothing before colon (like :/path) or alphanumeric word before colon
        # and there's content after the colon, it's likely a remote path
        if before_colon and "/" in path[colon_idx:] and len(path[colon_idx:]) > 1:
            return True

    return False


async def _prepare_bam_files_async(
    bams: List[str],
    ssh_password: Optional[str] = None,
    ssh_server: Optional[str] = None
) -> List[str]:
    """Async version: fetches remote files via thread pool, leaves local paths untouched."""
    result = []
    fm = FileManager()

    for bam_path in bams:
        # If ssh_server is provided, construct the full SCP path
        if ssh_server and not _is_remote_path(bam_path):
            remote_path = f"{ssh_server}:{bam_path}"
            if not ssh_password:
                raise ValueError(f"SSH password required for remote file: {remote_path}")
            print_log(f"Constructing SCP path: {remote_path}")
            local_path = await asyncio.get_event_loop().run_in_executor(
                get_cli_executor(),
                lambda p=remote_path, sp=ssh_password: fm.fetch_remote_file(p, sp)
            )
            print_log(f"Successfully fetched remote file to: {local_path}")
            result.append(local_path)
        elif _is_remote_path(bam_path):
            print_log(f"Detected as remote path, fetching via SCP: {bam_path}")
            if not ssh_password:
                raise ValueError(f"SSH password required for remote file: {bam_path}")
            local_path = await asyncio.get_event_loop().run_in_executor(
                get_cli_executor(),
                lambda p=bam_path, sp=ssh_password: fm.fetch_remote_file(p, sp)
            )
            print_log(f"Successfully fetched remote file to: {local_path}")
            result.append(local_path)
        else:
            print_log(f"Local path, using as-is: {bam_path}")
            result.append(bam_path)

    print_log(f"Final prepared BAM files: {result}")
    return result


def _prepare_bam_files(
    bams: List[str],
    ssh_password: Optional[str] = None,
    ssh_server: Optional[str] = None
) -> List[str]:
    """
    Prepare BAM files for processing.

    If ssh_server is provided, the BAM paths are treated as remote paths on that server.
    If file path is a remote SCP path (user@host:/path), fetch it via SSH/SCP.
    Otherwise, it's assumed to be a local file.

    Args:
        bams: List of BAM file paths (local or remote)
        ssh_password: SSH password for remote files
        ssh_server: SSH server address (user@host) for remote files

    Returns:
        List of local file paths
    """
    print_log(f"Preparing BAM files: {bams}")
    print_log(f"SSH password provided: {'Yes' if ssh_password else 'No'}")
    print_log(f"SSH server provided: {ssh_server}")
    result = []
    file_manager = FileManager()

    for bam_path in bams:
        print_log(f"Checking BAM path: {bam_path}")

        # If ssh_server is provided, construct the full SCP path
        if ssh_server and not _is_remote_path(bam_path):
            # Convert local-style path to SCP format
            remote_path = f"{ssh_server}:{bam_path}"
            print_log(f"Constructing SCP path from ssh_server: {remote_path}")
            if not ssh_password:
                raise ValueError(
                    f"SSH password required for remote file: {remote_path}")
            local_path = file_manager.fetch_remote_file(
                remote_path=remote_path,
                ssh_password=ssh_password,
            )
            print_log(f"Successfully fetched remote file to: {local_path}")
            result.append(local_path)
        elif _is_remote_path(bam_path):
            print_log(f"Detected as remote path, fetching via SCP: {bam_path}")
            # Remote file in user@host:/path format - fetch via SCP
            if not ssh_password:
                raise ValueError(
                    f"SSH password required for remote file: {bam_path}")
            local_path = file_manager.fetch_remote_file(
                remote_path=bam_path,
                ssh_password=ssh_password,
            )
            print_log(f"Successfully fetched remote file to: {local_path}")
            result.append(local_path)
        else:
            print_log(f"Local path, using as-is: {bam_path}")
            # Local file - use as-is
            result.append(bam_path)

    print_log(f"Final prepared BAM files: {result}")
    return result


@api_router.post(
    "/{tool_name}/execute",
    response_model=APIResponse,
    summary="Execute a CLI tool",
    description="Executes a specified CLI tool with the given arguments and returns the result",
)
async def execute_tool(
    tool_name: str,
    request: ToolExecutionRequest,
) -> APIResponse:
    """Execute a CLI tool with provided arguments"""
    print_log(f"=== routers.py - execute_tool (POST) START ===")
    print_log(f"HTTP Method: POST")
    print_log(f"=== Received execution request for tool: {tool_name} ===")
    print_log(f"Raw request object: {type(request)}")
    print_log(f"Raw request dict: {json.dumps(request.__dict__, indent=2)}")

    # Find the module path for the tool
    module_path = None
    for name, module in settings.TOOLS_CONFIG:
        if name == tool_name:
            module_path = module
            break

    # Check binary tools as fallback
    if not module_path and tool_name in settings.TOOL_BINARIES:
        module_path = settings.TOOL_BINARIES[tool_name]

    if not module_path:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Tool '{tool_name}' not found",
        )

    # Prepare BAM files (fetch remote if needed)
    # Handle different request formats:
    # - Format 1 (current): { tool_name: '...', bams: [...], channel_tag: 'ch', min_rq: 0.8, ssh_server: '...', ssh_password: '...' }
    # - Format 2 (old): { tool_name: '...', args: { bams: [...], channel_tag: 'ch', ... } }
    print_log(f"=== routers.py - execute_tool() START ===")
    print_log(f"Tool name: {tool_name}")
    print_log(f"Request type: {type(request).__name__}")

    # Determine the request format and extract arguments
    # We need to detect whether the request used args dict format or top-level format
    # by checking if non-tool_name fields were provided at the top level

    # Get all fields that are not tool_name and not in the base model
    # Note: In Pydantic v2, extra fields are stored in __pydantic_extra__, not __dict__
    non_default_fields = {
        k: v for k, v in request.__dict__.items() if k not in ('tool_name', 'args')}

    # Check __pydantic_extra__ for top-level format (Pydantic v2 extra fields)
    if hasattr(request, '__pydantic_extra__') and request.__pydantic_extra__:
        # Top-level format: extra fields were provided directly
        print_log(
            f"Request format: top-level fields format (via __pydantic_extra__)")
        print_log(
            f"Top-level fields: {json.dumps(request.__pydantic_extra__, indent=2)}")
        args_dict = dict(request.__pydantic_extra__)
    elif non_default_fields:
        # Fallback for Pydantic v1 behavior or if extra fields somehow in __dict__
        print_log(f"Request format: top-level fields format")
        print_log(
            f"Top-level fields: {json.dumps(non_default_fields, indent=2)}")
        args_dict = dict(non_default_fields)
    elif hasattr(request, 'args') and request.args:
        # Args dict format: args contains data
        print_log(f"Request format: args dict format")
        print_log(f"request.args value: {json.dumps(request.args, indent=2)}")
        args_dict = dict(request.args)
    else:
        # No arguments provided
        print_log(f"Request format: no arguments")
        args_dict = {}

    print_log(f"=== Extracted args ===")
    print_log(f"args_dict: {json.dumps(args_dict, indent=2)}")

    # Extract ssh_password and ssh_server from args if present
    ssh_password = args_dict.pop("ssh_password", None)
    ssh_server = args_dict.pop("ssh_server", None)
    print_log(f"SSH password: {'SET' if ssh_password else 'NOT SET'}")
    print_log(f"SSH server: {ssh_server or 'NOT SET'}")

    # Prepare file-type arguments (fetch remote if needed)
    # Dynamically detect file-type args from ARGUMENT_SCHEMAS
    file_arg_keys = []
    if tool_name and tool_name in ARGUMENT_SCHEMAS:
        for arg in ARGUMENT_SCHEMAS[tool_name].get("arguments", []):
            if arg.get("type") == "file":
                file_arg_keys.append(arg["name"])

    # Also handle legacy "bams" key for existing Python tools
    if "bams" in args_dict and "bams" not in file_arg_keys:
        file_arg_keys.append("bams")

    if file_arg_keys:
        for key in file_arg_keys:
            if key in args_dict:
                val = args_dict[key]
                if not isinstance(val, list):
                    val = [val]
                try:
                    prepared = await _prepare_bam_files_async(val, ssh_password, ssh_server)
                    args_dict[key] = prepared[0] if len(prepared) == 1 and not isinstance(args_dict.get(key), list) else prepared
                    print_log(f"Prepared {key}: {args_dict[key]}")
                except Exception as e:
                    print_log(f"File preparation failed for {key}: {str(e)}")
                    return APIResponse(
                        success=False,
                        error="File preparation failed",
                        message=str(e),
                    )

    # Clean up old file outputs
    CLIRunner.cleanup_file_outputs()

    # Execute the CLI module (non-blocking via thread pool)
    print_log(f"Executing CLI with args: {json.dumps(args_dict, indent=2)}")

    # For binary tools, dispatch via run_binary instead of run_cli_module
    if tool_name in settings.TOOL_BINARIES:
        bin_path = settings.TOOL_BINARIES[tool_name]
        args_list = CLIRunner._dict_to_args(tool_name, args_dict)
        print_log(f"Executing binary: {' '.join([bin_path] + args_list)}")
        returncode, stdout, stderr, command = await CLIRunner.run_binary_async(
            bin_path, args_list, timeout=settings.CLI_TOOL_TIMEOUT,
        )
    else:
        returncode, stdout, stderr, command = await CLIRunner.run_cli_with_json_async(
            module_path=module_path,
            json_args=json.dumps(args_dict),
            timeout=settings.CLI_TOOL_TIMEOUT,
        )
    print_log(f"CLI command executed: {' '.join(command)}")
    print_log(f"CLI exit code: {returncode}")
    if stdout:
        print_log(f"CLI stdout:\n{stdout}")
    if stderr:
        print_log(f"CLI stderr:\n{stderr}")

    # Clean up temporary files
    FileManager.cleanup_temp_files()

    # Auto-register output files for asrtc tool
    if tool_name == "asrtc" and args_dict.get("prefix"):
        prefix = args_dict["prefix"]
        if isinstance(prefix, list):
            prefix = prefix[0] if prefix else ""
        if prefix:
            output_file = f"{prefix}.asrtc.txt"
            if os.path.exists(output_file):
                FileRegistry.register([{"name": output_file, "path": output_file}])
                print_log(f"Auto-registered output file: {output_file}")

    # Get file outputs if available
    file_outputs = CLIRunner.get_registered_files()

    # Build response
    response = ToolExecutionResponse(
        success=(returncode == 0),
        tool_name=tool_name,
        stdout=stdout,
        stderr=stderr,
        exit_code=returncode,
        command=command if command else None,
        file_outputs=file_outputs if file_outputs else None,
    )

    if returncode == 0:
        return APIResponse(success=True, data=response.dict())
    else:
        return APIResponse(
            success=False,
            error="Execution failed",
            message=f"Tool '{tool_name}' exited with code {returncode}",
            data={
                "stdout": stdout,
                "stderr": stderr,
                "command": " ".join(command) if command else None,
            } if stdout or stderr else None,
        )


@api_router.get(
    "/{tool_name}/execute",
    response_model=APIResponse,
    summary="Execute a CLI tool (GET)",
    description="Executes a specified CLI tool using GET parameters",
)
async def execute_tool_get(
    tool_name: str,
    bams: List[str] = None,
    channel_tag: str = None,
    min_rq: float = None,
    ssh_password: Optional[str] = None,
    ssh_server: Optional[str] = None,
) -> APIResponse:
    """Execute a CLI tool with query parameters"""
    print_log(f"=== routers.py - execute_tool_get (GET) START ===")
    print_log(f"HTTP Method: GET")
    print_log(f"Tool name: {tool_name}")
    print_log(
        f"Query params: bams={bams}, channel_tag={channel_tag}, min_rq={min_rq}, ssh_server={ssh_server}")
    # Find the module path for the tool
    module_path = None
    for name, module in settings.TOOLS_CONFIG:
        if name == tool_name:
            module_path = module
            break

    if not module_path:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Tool '{tool_name}' not found",
        )

    # Clean up old file outputs
    CLIRunner.cleanup_file_outputs()

    # Build arguments from query parameters (non-blocking)
    args_dict = {}
    if bams:
        args_dict["bams"] = await _prepare_bam_files_async(bams, ssh_password, ssh_server)
    if channel_tag:
        args_dict["channel_tag"] = channel_tag
    if min_rq is not None:
        args_dict["min_rq"] = min_rq

    # Execute the CLI module (non-blocking via thread pool)
    returncode, stdout, stderr, command = await CLIRunner.run_cli_with_json_async(
        module_path=module_path,
        json_args=json.dumps(args_dict),
        timeout=settings.CLI_TOOL_TIMEOUT,
    )

    # Clean up temporary files
    FileManager.cleanup_temp_files()

    # Get file outputs if available
    file_outputs = CLIRunner.get_registered_files()

    # Build response
    response = ToolExecutionResponse(
        success=(returncode == 0),
        tool_name=tool_name,
        stdout=stdout,
        stderr=stderr,
        exit_code=returncode,
        command=command if command else None,
        file_outputs=file_outputs if file_outputs else None,
    )

    if returncode == 0:
        return APIResponse(success=True, data=response.dict())
    else:
        return APIResponse(
            success=False,
            error="Execution failed",
            message=f"Tool '{tool_name}' exited with code {returncode}",
            data={
                "stdout": stdout,
                "stderr": stderr,
                "command": " ".join(command) if command else None,
            } if stdout or stderr else None,
        )
