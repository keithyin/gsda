"""API Routers - Define all API endpoints for CLI tools"""

import json
import os
import sys
import tempfile
import uuid
from typing import List, Optional
from fastapi import APIRouter, HTTPException, status, UploadFile, File
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from gseda.server.config import settings
from gseda.server.core.runners import CLIRunner, ARGUMENT_SCHEMAS
from gseda.server.core.schema import (
    ToolExecutionRequest,
    ToolExecutionResponse,
    ToolListResponse,
    ToolInfoResponse,
    APIResponse,
)
from gseda.server.core.files import FileManager


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
                raise ValueError(f"SSH password required for remote file: {remote_path}")
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
                raise ValueError(f"SSH password required for remote file: {bam_path}")
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
    non_default_fields = {k: v for k, v in request.__dict__.items() if k not in ('tool_name', 'args')}

    # Check __pydantic_extra__ for top-level format (Pydantic v2 extra fields)
    if hasattr(request, '__pydantic_extra__') and request.__pydantic_extra__:
        # Top-level format: extra fields were provided directly
        print_log(f"Request format: top-level fields format (via __pydantic_extra__)")
        print_log(f"Top-level fields: {json.dumps(request.__pydantic_extra__, indent=2)}")
        args_dict = dict(request.__pydantic_extra__)
    elif non_default_fields:
        # Fallback for Pydantic v1 behavior or if extra fields somehow in __dict__
        print_log(f"Request format: top-level fields format")
        print_log(f"Top-level fields: {json.dumps(non_default_fields, indent=2)}")
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

    if "bams" in args_dict:
        print_log(f"Processing BAM files: {args_dict['bams']}")
        try:
            args_dict["bams"] = _prepare_bam_files(args_dict["bams"], ssh_password, ssh_server)
            print_log(f"Local BAM files after SCP fetch: {args_dict['bams']}")
        except Exception as e:
            print_log(f"BAM preparation failed: {str(e)}")
            return APIResponse(
                success=False,
                error="BAM preparation failed",
                message=str(e),
            )

    # Execute the CLI module
    print_log(f"Executing CLI with args: {json.dumps(args_dict, indent=2)}")
    returncode, stdout, stderr, command = CLIRunner.run_cli_with_json(
        module_path=module_path,
        json_args=json.dumps(args_dict),
    )
    print_log(f"CLI command executed: {' '.join(command)}")
    print_log(f"CLI exit code: {returncode}")
    if stdout:
        print_log(f"CLI stdout:\n{stdout}")
    if stderr:
        print_log(f"CLI stderr:\n{stderr}")

    # Clean up temporary files
    FileManager.cleanup_temp_files()

    # Build response
    response = ToolExecutionResponse(
        success=(returncode == 0),
        tool_name=tool_name,
        stdout=stdout,
        stderr=stderr,
        exit_code=returncode,
        command=command if command else None,
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
    print_log(f"Query params: bams={bams}, channel_tag={channel_tag}, min_rq={min_rq}, ssh_server={ssh_server}")
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

    # Build arguments from query parameters
    args_dict = {}
    if bams:
        args_dict["bams"] = _prepare_bam_files(bams, ssh_password, ssh_server)
    if channel_tag:
        args_dict["channel_tag"] = channel_tag
    if min_rq is not None:
        args_dict["min_rq"] = min_rq

    # Execute the CLI module
    returncode, stdout, stderr, command = CLIRunner.run_cli_with_json(
        module_path=module_path,
        json_args=json.dumps(args_dict),
    )

    # Clean up temporary files
    FileManager.cleanup_temp_files()

    # Build response
    response = ToolExecutionResponse(
        success=(returncode == 0),
        tool_name=tool_name,
        stdout=stdout,
        stderr=stderr,
        exit_code=returncode,
        command=command if command else None,
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
