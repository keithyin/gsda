#!/usr/bin/env python3
"""
Verification script for GSEDA Server

This script checks if all required files and dependencies are in place.
"""

import os
from pathlib import Path


def check_file_exists(path: str, description: str) -> bool:
    """Check if a file exists and report status."""
    exists = os.path.exists(path)
    status = "✓" if exists else "✗"
    print(f"{status} {description}: {path}")
    return exists


def main():
    """Run verification checks."""
    server_dir = Path(__file__).parent
    print("=" * 60)
    print("GSEDA Server Verification")
    print("=" * 60)
    print()

    # Check Python files
    print("Checking Python source files:")
    python_files = [
        ("config.py", "Configuration settings"),
        ("main.py", "FastAPI application entry point"),
        ("api/routers.py", "API route definitions"),
        ("core/runners.py", "CLI execution runner"),
        ("core/schema.py", "Pydantic schemas"),
    ]

    all_present = True
    for filename, desc in python_files:
        path = server_dir / filename
        if not check_file_exists(str(path), desc):
            all_present = False
    print()

    # Check templates
    print("Checking templates:")
    html_file = server_dir / "templates" / "index.html"
    if check_file_exists(str(html_file), "Main HTML template"):
        file_size = os.path.getsize(html_file)
        print(f"  File size: {file_size:,} bytes")
    print()

    # Check frontend files
    print("Checking Vue 3 frontend:")
    frontend_files = [
        ("frontend/src/main.ts", "Vue entry point"),
        ("frontend/src/App.vue", "Main application component"),
        ("frontend/src/components/ToolForm.vue", "Tool form component"),
        ("frontend/src/components/ResultViewer.vue", "Results viewer component"),
        ("frontend/src/api.ts", "API client"),
    ]

    for filename, desc in frontend_files:
        path = server_dir / filename
        if not check_file_exists(str(path), desc):
            all_present = False

    # Check build config
    print()
    print("Checking build configuration:")
    check_file_exists(str(server_dir / "frontend" / "package.json"), "npm package.json")
    check_file_exists(str(server_dir / "frontend" / "vite.config.ts"), "Vite configuration")
    print()

    # Check requirements
    print("Checking dependencies:")
    req_file = server_dir / "requirements.txt"
    if os.path.exists(req_file):
        size = os.path.getsize(req_file)
        print(f"✓ requirements.txt: {size:,} bytes")
    else:
        print(f"✗ requirements.txt: Missing")
    print()

    # Summary
    print("=" * 60)
    if all_present:
        print("✓ All required files are present!")
    else:
        print("✗ Some files are missing. Please review the output above.")
    print("=" * 60)


if __name__ == "__main__":
    main()
