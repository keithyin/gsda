# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the `gsda` repository containing the `gseda` Python package for genomic data analysis CLI tools, plus a FastAPI-based web service that exposes these tools via a REST API and web UI.

### High-Level Architecture

```
gsda/
├── src/gseda/                    # Main Python package with CLI tools
│   ├── bam_ana/                  # BAM analysis tools
│   ├── bam_filter/               # BAM filtering tools
│   ├── bam_surgery/              # BAM manipulation tools
│   ├── fastx_ana/                # FASTX analysis tools
│   ├── msa_view/                 # Multiple sequence alignment viewer
│   ├── ppl/                      # Population genetics analysis tools
│   └── ...
└── src/gseda/server/             # Web service (FastAPI + Vue 3)
    ├── main.py                   # FastAPI application entry point
    ├── config.py                 # Configuration (lists all CLI tools)
    ├── api/routers.py            # REST API endpoints
    ├── core/                     # Core functionality
    │   ├── runners.py            # CLI module executors
    │   ├── schema.py             # Pydantic validation models
    │   └── files.py              # File management (SCP, uploads)
    ├── frontend/                 # Vue 3 source
    └── templates/                # HTML templates (CDN fallback)
```

## Development Commands

### Backend (FastAPI Server)

```bash
# Install server dependencies
cd src/gseda/server
pip install -r requirements.txt

# Run server (development mode with CDN-based frontend)
uvicorn gseda.server.main:app --reload --host 0.0.0.0 --port 8000

# Or use the startup script
python gseda.server.start_server:main

# Run server via console entry point
gseda-server
```

### Frontend (Vue 3)

```bash
# Install frontend dependencies
cd src/gseda/server/frontend
npm install

# Development server
npm run dev

# Production build
npm run build
```

### Testing

```bash
# Run tests (from repository root)
pytest tests/

# Run specific test file
pytest tests/test_api_routers.py -v

# Run with coverage
pytest tests/ --cov=gseda
```

### Linting (if configured)

```bash
# Check for any linting or formatting tools
# (No linting tools currently configured in this project)
```

## API Documentation

Interactive API docs are available at:
- http://localhost:8000/docs (Swagger UI)
- http://localhost:8000/redoc (ReDoc)

### Key Endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/api/tools/list` | List all available CLI tools |
| GET | `/api/tools/{name}/info` | Get tool metadata and argument schema |
| POST | `/api/tools/{name}/execute` | Execute a CLI tool (JSON body) |
| GET | `/api/tools/{name}/execute` | Execute a CLI tool (query params) |

### Example API Call

```bash
curl -X POST "http://localhost:8000/api/tools/bam-basic-stat/execute" \
  -H "Content-Type: application/json" \
  -d '{
    "tool_name": "bam-basic-stat",
    "args": {
      "bams": ["/path/to/file.bam"],
      "channel_tag": "ch",
      "min_rq": 0.8
    }
  }'
```

## Adding New CLI Tools to the Web Service

1. **Edit `src/gseda/server/config.py`**: Add the tool to `TOOLS_CONFIG`:
   ```python
   ("tool-name", "gseda.module.path:main_cli"),
   ```

2. **Optional - Add form schema**: In `core/schema.py`, define a request model if needed:
   ```python
   class NewToolRequest(ToolRequestBase):
       param1: str
       param2: int
   ```

3. **Optional - Add form component**: Create `src/gseda/server/frontend/src/components/forms/NewToolForm.vue` and register it in `ToolForm.vue`.

## Key Design Patterns

### CLI Tool Execution Flow

1. Frontend sends request to `/api/tools/{name}/execute`
2. `routers.py` receives the request and extracts arguments
3. `_prepare_bam_files()` handles remote file fetching via SCP if needed
4. `CLIRunner.run_cli_with_json()` converts JSON args to CLI format
5. `subprocess.run()` executes the CLI module
6. Output is returned as JSON response

### SSH/SCP Remote File Handling

- Remote paths are detected by pattern matching (e.g., `host:/path` or `user@host:/path`)
- `FileManager.fetch_remote_file()` uses paramiko to download files via SFTP
- Files are stored in system temp directory and cleaned up after execution

### Vue 3 Frontend Architecture

- `App.vue` - Main layout with sidebar (tool list) and content area
- `ToolForm.vue` - Container that dynamically renders tool-specific forms
- `ResultViewer.vue` - Displays execution results with copy-to-clipboard
- Forms for specific tools live in `components/forms/` (e.g., `BAMBasicStatForm.vue`)

## Configuration

Environment variables for `gseda/server`:

| Variable | Default | Description |
|----------|---------|-------------|
| `DEBUG` | `false` | Enable debug mode |
| `CLI_TOOL_TIMEOUT` | `3600` | Timeout for CLI execution in seconds |

## Known Limitations

1. **File Uploads**: Uses file paths as parameters; no file upload UI yet
2. **Long-Running Tasks**: No task queue; long tasks may timeout
3. **Authentication**: No user authentication or permissions
4. **Form Coverage**: Only some tools have dedicated frontend forms; all are API-accessible

## Related Documentation

- `src/gseda/server/README.md` - Full server documentation
- `src/gseda/server/QUICKSTART.md` - Quick start guide
- `src/gseda/server/IMPLEMENTATION_SUMMARY.md` - Implementation details
- `src/gseda/server/design.md` - Design document

## Project Metadata

- **Package**: `gseda` (version 1.21.10)
- **CLI tools**: 22 tools covering BAM analysis, sequencing reports, quality stats, etc.
- **Backend**: FastAPI + Uvicorn
- **Frontend**: Vue 3 + Element Plus + Vite
- **Python**: >=3.8


# claude 权限

1. cluade 在执行的过程中 允许所有的文件编辑操作，无需请求授权

# cluade 说明

Although you are a LLM, you can't perform computer operations by yourself. However, you are powering Claude code, and you can use its tools and skills to operate the current computer and external entities.

# cluade 运行说明

cluade 每次运行前，要阅读