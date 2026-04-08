#!/usr/bin/env python3
"""
Startup script for GSEDA Server

This script provides a simple way to start the web server.
"""

import subprocess
import sys
from pathlib import Path
server_dir = Path(__file__).parent
project_root = server_dir.parent.parent

# Add src directory to Python path for gseda package imports
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))


def main():
    """Start the GSEDA web server."""
    # Determine paths and add to sys.path for module imports

    # Change to server directory
    import os
    os.chdir(server_dir)

    # Check if frontend is built
    frontend_dist = server_dir / "static" / "dist"
    if not frontend_dist.exists():
        print("Warning: Frontend not built. Please run 'npm run build' in the frontend directory.")
        print("Proceeding with backend-only mode...")

    # Start uvicorn
    cmd = [
        sys.executable,
        "-m",
        "uvicorn",
        "gseda.server.main:app",
        "--reload" if sys.argv[1:] and "--debug" in sys.argv else "",
        "--host",
        "0.0.0.0",
        "--port",
        "8000"
    ]

    print(f"Starting GSEDA Server on http://0.0.0.0:8000")
    print("Press CTRL+C to stop\n")

    subprocess.run(cmd)


if __name__ == "__main__":
    main()
