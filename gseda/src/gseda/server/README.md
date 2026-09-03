# GSEDA Server - Web Service for CLI Tools

Web-based interface for genomic data analysis CLI tools.

## Overview

This project provides a web service that exposes the CLI tools from the gseda package through a modern REST API and an elegant web interface.

## Features

- **RESTful API**: All CLI tools are accessible via HTTP endpoints
- **Web Interface**: User-friendly UI for executing tools and viewing results
- **FastAPI Backend**: High-performance asynchronous server
- **Vue 3 Frontend**: Modern, reactive user interface
- **Extensible Architecture**: Easy to add new CLI tools

## Project Structure

```
gseda/src/gseda/server/
├── __init__.py          # Package marker
├── main.py              # FastAPI application entry point
├── config.py            # Configuration settings
├── requirements.txt     # Python dependencies
├── README.md            # This file
├── api/                 # API routes and endpoints
│   ├── __init__.py
│   └── routers.py       # Tool execution endpoints
├── core/                # Core functionality
│   ├── __init__.py
│   ├── runners.py       # CLI module executors
│   └── schema.py        # Pydantic models
├── templates/           # HTML templates
│   └── index.html       # Main web interface (basic)
├── static/              # Static assets
│   └── dist/            # Compiled Vue frontend (after build)
└── frontend/            # Vue 3 source code
    ├── package.json
    ├── vite.config.js
    ├── tsconfig.json
    └── src/
        ├── main.ts
        ├── App.vue
        ├── router.ts
        ├── api.ts
        ├── types.ts
        ├── components/
        │   ├── ToolForm.vue
        │   ├── ResultViewer.vue
        │   └── forms/
        │       ├── BAMBasicStatForm.vue
        │       └── MSAViewForm.vue
        └── router.ts
```

## Installation

### Dependencies

Install Python and Node.js dependencies:

```bash
cd gseda/src/gseda/server
pip install -r requirements.txt

cd frontend
npm install
```

### Building the Frontend

Build the Vue 3 frontend:

```bash
cd gseda/src/gseda/server/frontend
npm run build
```

The build output will be placed in `../static/dist/`.

## Running the Server

### Development Mode (Backend Only)

Start the FastAPI server (uses CDN-based frontend):

```bash
cd gseda/src/gseda/server
uvicorn gseda.server.main:app --reload --host 0.0.0.0 --port 8000
```

Then visit `http://localhost:8000` in your browser.

**Note**: This mode uses Vue.js from CDN and is intended for quick development/testing.
For production, build the frontend with `npm run build`.

### Production Mode (Built Frontend)

1. Build the frontend:

```bash
cd gseda/src/gseda/server/frontend
npm run build
```

2. Start the backend:

```bash
cd gseda/src/gseda/server
uvicorn gseda.server.main:app --host 0.0.0.0 --port 8000
```

3. Visit `http://localhost:8000`

## API Documentation

Interactive API documentation is available at `/docs`:

```bash
# Start server and visit
uvicorn gseda.server.main:app --reload
# Then open: http://localhost:8000/docs
```

### Available Endpoints

- `GET /api/tools/list` - List all available tools
- `GET /api/tools/{tool_name}/info` - Get tool metadata
- `POST /api/tools/{tool_name}/execute` - Execute a CLI tool (JSON body)
- `GET /api/tools/{tool_name}/execute` - Execute a CLI tool (query params)

### Example API Usage

```bash
# List all tools
curl http://localhost:8000/api/tools/list

# Get info for bam-basic-stat
curl http://localhost:8000/api/tools/bam-basic-stat/info

# Execute bam-basic-stat with JSON body
curl -X POST "http://localhost:8000/api/tools/bam-basic-stat/execute" \
  -H "Content-Type: application/json" \
  -d '{"tool_name": "bam-basic-stat", "args": {"bams": ["/path/to/file.bam"], "channel_tag": "ch"}}'

# Execute with query parameters
curl "http://localhost:8000/api/tools/bam-basic-stat/execute?bams=/path/to/file.bam&channel_tag=ch"
```

## Adding New CLI Tools

To add a new CLI tool to the web interface:

1. **Register in `config.py`**: Add the tool to `TOOLS_CONFIG`:

```python
TOOLS_CONFIG = [
    # ... existing tools
    ("new-tool-name", "gseda.module.path:main_cli"),
]
```

2. **Create Form Component** (optional): Add a form component in `frontend/src/components/forms/`

3. **Add to Router**: Update `ToolForm.vue` to include the new form

4. **Test**: Build frontend and test via web interface

## Configuration

Environment variables:

- `DEBUG`: Enable debug mode (`true` or `false`)
- `CLI_TOOL_TIMEOUT`: Timeout for CLI execution in seconds (default: 3600)

## Technologies Used

### Backend
- FastAPI - Modern web framework
- Pydantic - Data validation
- Jinja2 - Template engine
- Uvicorn - ASGI server

### Frontend
- Vue 3 - Progressive JavaScript framework
- Vite - Next generation frontend tooling
- Element Plus - Vue 3 UI library
- Axios - HTTP client

## License

MIT License - See the parent project for details.
