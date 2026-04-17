"""
GSEDA Server - FastAPI Web Service for CLI Tools

This module provides a web interface for executing genomic data analysis CLI tools.
"""

from gseda.server.api import routers
from gseda.server.config import settings
import os
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from jinja2 import Environment, FileSystemLoader
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
import sys

from pathlib import Path
server_dir = Path(__file__).parent
project_root = server_dir.parent.parent

# Add src directory to Python path for gseda package imports
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))


# Create FastAPI application
app = FastAPI(
    title=settings.APP_NAME,
    description="""
## GSEDA CLI Tools Web Service

This service provides web-based access to genomic data analysis command-line tools.

### Features

- **Web Interface**: Access all CLI tools through a modern web UI
- **REST API**: Execute tools via HTTP requests
- **Real-time Output**: See results as they are generated
- **File Uploads**: Support for uploading input files

### Available Tools

The following CLI tools are available:
{tools}

### Usage

See the `/docs` endpoint for interactive API documentation, or use the web interface.
    """.format(
        tools="\n".join([f"- `{name}`" for name, _ in settings.TOOLS_CONFIG])
    ),
    version=settings.APP_VERSION,
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include API routers
app.include_router(routers.api_router, prefix=settings.API_PREFIX)

# Configure templates (Jinja2 with Vue.js-compatible delimiters)
templates_dir = settings.TEMPLATE_DIR
os.makedirs(templates_dir, exist_ok=True)

env = Environment(
    loader=FileSystemLoader(templates_dir),
    autoescape=False,
    variable_start_string='{#',
    variable_end_string='#}',
    block_start_string='{#%',
    block_end_string='%#}',
    comment_start_string='{##',
    comment_end_string='##}',
)

templates = Jinja2Templates(env=env)

# Mount static files (for Vue frontend after build)
static_dir = settings.STATIC_DIR
os.makedirs(static_dir, exist_ok=True)

# Create custom StaticFiles that serves index.html for SPA routing
class SPAStaticFiles(StaticFiles):
    """StaticFiles that serves index.html for unknown paths (SPA fallback)"""

    async def get_response(self, path: str, scope):
        try:
            return await super().get_response(path, scope)
        except Exception:
            # If file not found, serve index.html for SPA routing
            return await self.get_response("index.html", scope)

if os.path.exists(settings.FRONTEND_BUILD_DIR):
    app.mount(
        "/static", SPAStaticFiles(directory=settings.FRONTEND_BUILD_DIR, html=True), name="static")


@app.get("/")
async def root(request: Request):
    """Serve the main web application

   优先返回构建后的前端静态文件，如果不存在则返回模板文件作为回退
    """
    # Check if built frontend exists
    built_index = os.path.join(settings.FRONTEND_BUILD_DIR, "index.html")
    if os.path.exists(built_index):
        # Serve built frontend (will be handled by static files middleware)
        return FileResponse(built_index)

    # Fallback to template if not built
    return templates.TemplateResponse(
        "index.html", {
            "request": request,
            "title": settings.APP_NAME,
            # Initial placeholder values for Vue.js properties
            "selectedToolName": "",
            "selectedTool": None
        }
    )


@app.on_event("startup")
async def startup_event():
    """Run on server startup"""
    print(f"Starting {settings.APP_NAME} v{settings.APP_VERSION}")
    print(f"Available tools: {len(settings.TOOLS_CONFIG)}")
    for name, _ in settings.TOOLS_CONFIG:
        print(f"  - {name}")


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "version": settings.APP_VERSION}
