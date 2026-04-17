"""Pydantic schemas for API request/response validation"""

from typing import List, Optional, Any, Dict
from pydantic import BaseModel, Field


class ToolRequestBase(BaseModel):
    """Base schema for tool execution requests"""

    tool_name: str = Field(..., description="Name of the CLI tool to execute")


class BAMBasicStatRequest(ToolRequestBase):
    """Request schema for bam-basic-stat tool"""

    bams: List[str] = Field(..., description="Path(s) to BAM file(s)")
    channel_tag: Optional[str] = Field("ch", description="Channel tag: 'ch' or 'zm'")
    min_rq: Optional[float] = Field(None, description="Minimum read quality (0-1)")
    ssh_server: Optional[str] = Field(None, description="SSH server address for remote files")
    ssh_password: Optional[str] = Field(None, description="SSH password for remote files")


class MSAViewRequest(ToolRequestBase):
    """Request schema for msa-view tool"""

    bam: str = Field(..., description="Alignment BAM file path")
    ref_fasta: str = Field(..., description="Reference FASTA file path")
    ref_name: str = Field(..., description="Reference sequence name")
    start: Optional[int] = Field(None, description="Reference start position")
    end: Optional[int] = Field(None, description="Reference end position")
    o_fasta: Optional[str] = Field(None, description="Output FASTA file path")
    o_pic: Optional[str] = Field(None, description="Output image file path")


# Generic request that can be used for any tool
class ToolExecutionRequest(BaseModel):
    """Generic schema for executing any CLI tool"""

    tool_name: str = Field(..., description="Name of the CLI tool to execute")
    args: Dict[str, Any] = Field(default_factory=dict, description="Tool-specific arguments as JSON")

    class Config:
        """Pydantic config to allow extra fields from request body"""
        extra = "allow"


# Response schemas

class ToolExecutionResponse(BaseModel):
    """Response for successful tool execution"""

    success: bool
    tool_name: str
    stdout: str
    stderr: str
    exit_code: int
    result_data: Optional[Dict[str, Any]] = None
    command: Optional[List[str]] = Field(None, description="The executed command for debugging")


class ErrorResponse(BaseModel):
    """Response for failed tool execution"""

    success: bool = False
    error_type: str
    message: str
    tool_name: Optional[str] = None
    details: Optional[str] = None


class ToolListResponse(BaseModel):
    """Response listing all available tools"""

    tools: List[Dict[str, str]] = Field(..., description="List of available tools")


class ToolInfoResponse(BaseModel):
    """Response with tool metadata"""

    name: str
    module_path: str
    description: Optional[str] = None
    arguments: List[Dict[str, Any]] = Field(default_factory=list)


# Combined response model for API endpoints
class APIResponse(BaseModel):
    """Standard API response wrapper"""

    success: bool
    data: Optional[Any] = None
    error: Optional[str] = None
    message: Optional[str] = None
