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
│   ├── ppl/                      # Population genetics analysis tools
│   └── ...
└── src/gseda/server/             # Web service (FastAPI + Vue 3)
    ├── main.py                   # FastAPI application entry point
    ├── config.py                 # Configuration (lists all CLI tools)
    ├── api/routers.py            # REST API endpoints
    ├── core/                     # Core functionality
    │   ├── runners.py            # CLI module executors + FileRegistry
    │   ├── schema.py             # Pydantic validation models
    │   └── files.py              # File management (SCP, uploads)
    ├── frontend/                 # Vue 3 source
    └── templates/                # HTML templates (CDN fallback)
```

## Active CLI Tools (7)

Tools are defined in `src/gseda/server/config.py` `TOOLS_CONFIG`:

| Tool Name | Module Path |
|---|---|
| `reads-quality-stats-v3` | `gseda.ppl.reads_quality_stats_v3:main_cli` |
| `reads-quality-hp-tr` | `gseda.ppl.reads_quality_stats_hp_tr:main_cli` |
| `reads-quality-hp` | `gseda.ppl.reads_quality_stats_hp:main_cli` |
| `bam-basic-stat` | `gseda.bam_ana.bam_basic_stat:main_cli` |
| `low-q-analysis` | `gseda.ppl.low_q_analysis:main_cli` |
| `homo-and-str-ratio` | `gseda.ppl.homo_and_str_region_coverage:main_cli` |
| `macebell-ratio` | `gseda.ppl.macebell_ratio:main_cli` |

> `sequencing-report-v2` is currently commented out in config.

## Frontend Form Components

Located in `src/gseda/server/frontend/src/components/forms/`:

- `BAMBasicStatForm.vue`
- `LowQAnalysisForm.vue`
- `MacebellRatioForm.vue`
- `HomoAndStrRatioForm.vue`
- `ReadsQualityStatsV3Form.vue`
- `ReadsQualityStatsHPForm.vue`
- `ReadsQualityStatsHPTrForm.vue`
- `SequencingReportV2Form.vue`

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
| GET | `/api/tools/files/download/{filename}` | Download tool output file |
| POST | `/api/tools/files/remote-fetch` | Fetch remote file via SCP |
| POST | `/api/tools/files/upload` | Upload file from client |

### Example API Call

```bash
curl -X POST "http://localhost:8000/api/tools/bam-basic-stat/execute" \
  -H "Content-Type: application/json" \
  -d '{
    "tool_name": "bam-basic-stat",
    "bams": ["/path/to/file.bam"],
    "channel_tag": "ch",
    "min_rq": 0.8
  }'
```

## Adding New CLI Tools to the Web Service

1. **Edit `src/gseda/server/config.py`**: Add the tool to `TOOLS_CONFIG`:
   ```python
   ("tool-name", "gseda.module.path:main_cli"),
   ```

2. **Optional - Add argument schema**: In `core/runners.py` `ARGUMENT_SCHEMAS`, define the argument mapping:
   ```python
   "tool-name": {
       "arguments": [
           {"name": "bams", "type": "file", "required": True, "multiple": True, "positional": True},
           {"name": "min_rq", "type": "number", "required": False, "placeholder": "0.8"},
       ],
   },
   ```

3. **Optional - Add form component**: Create `src/gseda/server/frontend/src/components/forms/NewToolForm.vue` and register it in `ToolForm.vue`.

## Key Design Patterns

### CLI Tool Execution Flow

1. Frontend sends request to `/api/tools/{name}/execute`
2. `routers.py` receives the request and extracts arguments (supports both `args` dict and top-level field formats)
3. `_prepare_bam_files()` handles remote file fetching via SCP if needed
4. `CLIRunner.run_cli_with_json()` converts JSON args to CLI format using `ARGUMENT_SCHEMAS`
5. `subprocess.run()` executes the CLI module
6. `FileRegistry` tracks output files for download
7. Output (stdout, stderr, exit_code, file_outputs) is returned as JSON response

### File Registry (Output File Downloads)

- CLI tools register output files via `FileRegistry.register()` (signaled by `GSEDA_FILE_INDEX_PATH` env var)
- Files are indexed in a stable JSON file at `/tmp/gseda_current_file_index.json`
- Download endpoint `/api/tools/files/download/{filename}` serves registered files
- Entries are cleaned up after TTL (3600s) via `FileRegistry.cleanup()`

### SSH/SCP Remote File Handling

- Remote paths are detected by pattern matching (e.g., `host:/path` or `user@host:/path`)
- `FileManager.fetch_remote_file()` uses paramiko to download files via SFTP
- Files are stored in system temp directory and cleaned up after execution

### Vue 3 Frontend Architecture

- `App.vue` - Main layout with sidebar (tool list) and content area
- `ToolForm.vue` - Container that dynamically renders tool-specific forms (supports file source toggle: local, upload, SCP)
- `ResultViewer.vue` - Displays execution results with copy-to-clipboard and file download links
- Forms for specific tools live in `components/forms/`

## Configuration

Environment variables for `gseda/server`:

| Variable | Default | Description |
|---|---|---|
| `DEBUG` | `false` | Enable debug mode |
| `CLI_TOOL_TIMEOUT` | `3600` | Timeout for CLI execution in seconds |

## Known Limitations

1. **File Uploads**: Uses file paths as parameters; upload UI exists but limited integration
2. **Long-Running Tasks**: No task queue; long tasks may timeout
3. **Authentication**: No user authentication or permissions
4. **Form Coverage**: Only some tools have dedicated frontend forms; all are API-accessible

## Related Documentation

- `src/gseda/server/README.md` - Full server documentation
- `src/gseda/server/QUICKSTART.md` - Quick start guide
- `src/gseda/server/IMPLEMENTATION_SUMMARY.md` - Implementation details
- `src/gseda/server/design.md` - Design document

## Project Metadata

- **Package**: `gseda` (APP_VERSION: 1.0.0 in config.py)
- **CLI tools**: 7 active tools (BAM analysis, sequencing reports, quality stats, etc.)
- **Backend**: FastAPI + Uvicorn
- **Frontend**: Vue 3 + Element Plus + Vite
- **Python**: >=3.8
