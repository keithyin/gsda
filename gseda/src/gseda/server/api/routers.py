"""API Routers - Define all API endpoints for CLI tools"""

import json
from typing import List
from fastapi import APIRouter, HTTPException, status
from gseda.server.config import settings
from gseda.server.core.runners import CLIRunner, ARGUMENT_SCHEMAS
from gseda.server.core.schema import (
    ToolExecutionRequest,
    ToolExecutionResponse,
    ToolListResponse,
    ToolInfoResponse,
    APIResponse,
)


# Create router
api_router = APIRouter(prefix="/tools", tags=["tools"])


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


@api_router.post(
    "/{tool_name}/execute",
    response_model=APIResponse,
    summary="Execute a CLI tool",
    description="Executes a specified CLI tool with the given arguments and returns the result",
)
async def execute_tool(tool_name: str, request: ToolExecutionRequest) -> APIResponse:
    """Execute a CLI tool with provided arguments"""
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

    # Execute the CLI module
    returncode, stdout, stderr = CLIRunner.run_cli_with_json(
        module_path=module_path,
        json_args=json.dumps(request.args),
    )

    # Build response
    response = ToolExecutionResponse(
        success=(returncode == 0),
        tool_name=tool_name,
        stdout=stdout,
        stderr=stderr,
        exit_code=returncode,
    )

    if returncode == 0:
        return APIResponse(success=True, data=response.dict())
    else:
        return APIResponse(
            success=False,
            error="Execution failed",
            message=f"Tool '{tool_name}' exited with code {returncode}",
            data={"stdout": stdout, "stderr": stderr} if stdout or stderr else None,
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
) -> APIResponse:
    """Execute a CLI tool with query parameters"""
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
        args_dict["bams"] = bams
    if channel_tag:
        args_dict["channel_tag"] = channel_tag
    if min_rq is not None:
        args_dict["min_rq"] = min_rq

    # Execute the CLI module
    returncode, stdout, stderr = CLIRunner.run_cli_with_json(
        module_path=module_path,
        json_args=str(args_dict).replace("'", '"'),
    )

    # Build response
    response = ToolExecutionResponse(
        success=(returncode == 0),
        tool_name=tool_name,
        stdout=stdout,
        stderr=stderr,
        exit_code=returncode,
    )

    if returncode == 0:
        return APIResponse(success=True, data=response.dict())
    else:
        return APIResponse(
            success=False,
            error="Execution failed",
            message=f"Tool '{tool_name}' exited with code {returncode}",
            data={"stdout": stdout, "stderr": stderr} if stdout or stderr else None,
        )
