"""CLI Module Runners - Execute CLI commands via subprocess"""

import asyncio
import json
import logging
import os
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import List, Tuple, Any, Dict, Optional
from gseda.server.config import settings

# Thread pool for non-blocking CLI execution
_cli_executor: Optional[ThreadPoolExecutor] = None
_cli_executor_lock = threading.Lock()


def get_cli_executor() -> ThreadPoolExecutor:
    """Get or create the shared CLI executor (lazy init, thread-safe)."""
    global _cli_executor
    if _cli_executor is None:
        with _cli_executor_lock:
            if _cli_executor is None:
                _cli_executor = ThreadPoolExecutor(max_workers=8)
    return _cli_executor


def shutdown_cli_executor():
    """Shutdown the shared CLI executor on server stop."""
    global _cli_executor
    if _cli_executor is not None:
        _cli_executor.shutdown(wait=True)
        _cli_executor = None


# Configure logging for console output
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class FileRegistry:
    """Registry for CLI tool output files.

    Uses a JSON file as the index. The parent process creates the index file
    and signals its path to the child via GSEDA_FILE_INDEX_PATH env var.
    CLI tools call FileRegistry.register() to add entries.
    """

    TTL = 3600  # seconds
    STABLE_PATH = os.path.join(tempfile.gettempdir(), "gseda_current_file_index.json")

    @staticmethod
    def register(files: List[Dict[str, str]]) -> None:
        """Register output files. Called by CLI tools via env var signal."""
        path = os.environ.get("GSEDA_FILE_INDEX_PATH")
        if not path:
            return
        try:
            with open(path, "r") as f:
                index = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError):
            index = []
        now = time.time()
        for entry in files:
            if "name" not in entry:
                continue
            if "path" not in entry:
                entry["path"] = entry["name"]
            index.append({
                "name": entry["name"],
                "path": entry["path"],
                "timestamp": now,
                "pid": os.getpid(),
            })
        with open(path, "w") as f:
            json.dump(index, f, indent=2)

    @staticmethod
    def get_registered_files() -> List[Dict[str, str]]:
        """Return [{name, filename, download_url}] for existing files."""
        index = _load_file_index()
        results = []
        for entry in index:
            fpath = entry.get("path", "")
            if not os.path.exists(fpath):
                continue
            fname = os.path.basename(fpath)
            results.append({
                "name": entry.get("name", fname),
                "filename": fname,
                "download_url": f"/api/tools/files/download/{fname}",
            })
        return results

    @classmethod
    def cleanup(cls, ttl: Optional[int] = None) -> int:
        """Remove entries older than TTL from the stable index path."""
        if ttl is None:
            ttl = cls.TTL
        path = cls.STABLE_PATH
        if not os.path.exists(path):
            return 0
        try:
            with open(path, "r") as f:
                index = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError):
            return 0
        now = time.time()
        kept = [e for e in index if (now - e.get("timestamp", 0)) < ttl]
        removed = len(index) - len(kept)
        with open(path, "w") as f:
            json.dump(kept, f, indent=2)
        return removed

    @staticmethod
    def _save(path: str, data: list) -> None:
        """Atomically write data to path."""
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(data, f)
        os.replace(tmp, path)


def _load_file_index() -> List[Dict[str, str]]:
    """Load current file index entries."""
    path = FileRegistry.STABLE_PATH
    if not os.path.exists(path):
        return []
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        return []


class CLIRunner:
    """Execute CLI modules and handle their output"""

    @staticmethod
    def _build_env() -> dict:
        """Build environment dict for subprocess, including file index.

        Uses a stable path so the download endpoint can find it.
        """
        env = os.environ.copy()
        index_path = os.path.join(
            tempfile.gettempdir(), "gseda_current_file_index.json",
        )
        FileRegistry._save(index_path, [])
        env["GSEDA_FILE_INDEX_PATH"] = index_path
        return env

    @staticmethod
    def cleanup_file_outputs() -> int:
        """Clean up stale file index entries before new execution."""
        return FileRegistry.cleanup()

    @staticmethod
    def get_registered_files() -> List[Dict[str, str]]:
        """Get list of registered output files with download URLs."""
        return FileRegistry.get_registered_files()

    @staticmethod
    def run_cli_module(
        module_path: str, args: List[str], timeout: int = None
    ) -> Tuple[int, str, str, List[str]]:
        """
        Run a CLI module via subprocess.

        Args:
            module_path: Python module path (e.g., "gseda.bam_ana.bam_basic_stat:main_cli")
            args: Command line arguments to pass
            timeout: Execution timeout in seconds (default from settings)

        Returns:
            Tuple of (exit_code, stdout, stderr, command)
        """
        if timeout is None:
            timeout = settings.CLI_TOOL_TIMEOUT

        # Split module_path into module and callable
        if ":" in module_path:
            module_name, callable_name = module_path.rsplit(":", 1)
        else:
            module_name = module_path
            callable_name = "main_cli"

        # Build command
        cmd = [
            "python",
            "-c",
            f"import sys; from {module_name} import {callable_name}; {callable_name}()",
        ] + args

        try:
            env = CLIRunner._build_env()
            # Print to console for easy debugging
            print("\n" + "="*60, flush=True)
            print(f"[CLI] EXECUTING: {' '.join(cmd)}", flush=True)
            print("="*60 + "\n", flush=True)
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=settings.PROJECT_ROOT,
                env=env,
            )
            print(f"[CLI] EXIT CODE: {result.returncode}", flush=True)
            if result.stdout:
                print(f"[CLI] STDOUT:\n{result.stdout}", flush=True)
            if result.stderr:
                print(f"[CLI] STDERR:\n{result.stderr}", flush=True)
            print("="*60 + "\n", flush=True)
            return result.returncode, result.stdout, result.stderr, cmd
        except subprocess.TimeoutExpired:
            print(f"[CLI] TIMEOUT: Command timed out after {timeout} seconds", flush=True)
            return -1, "", f"Command timed out after {timeout} seconds", cmd
        except Exception as e:
            print(f"[CLI] EXCEPTION: {str(e)}", flush=True)
            return -1, "", str(e), cmd

    @staticmethod
    def parse_args_from_module(module_path: str) -> Dict[str, Any]:
        """
        Parse argument definitions from a CLI module dynamically.
        Returns a dict of argument specifications.

        Note: This is a basic implementation. More sophisticated parsing
        would require importing the module and inspecting its parser.
        """
        # For now, return empty dict - UI will use predefined schemas
        # In future version, we could introspect argparse parsers
        return {}

    @staticmethod
    def run_cli_with_json(
        module_path: str, json_args: str, timeout: int = None
    ) -> Tuple[int, str, str, List[str]]:
        """
        Run a CLI module with JSON-formatted arguments.

        Args:
            module_path: Python module path
            json_args: JSON string containing argument key-value pairs
            timeout: Execution timeout in seconds

        Returns:
            Tuple of (exit_code, stdout, stderr, command)
        """
        print(f"[CLI RUNNER] Parsed args from JSON: {json_args}")

        # Determine tool name from module path for correct argument mapping
        tool_name = None
        for name, path in settings.TOOLS_CONFIG:
            if path == module_path:
                tool_name = name
                break

        args_dict = json.loads(json_args)
        args_list = CLIRunner._dict_to_args(tool_name, args_dict)
        print(f"[CLI RUNNER] Converted to CLI args: {args_list}")

        # Dispatch to binary if this is a binary tool
        if tool_name and tool_name in settings.TOOL_BINARIES:
            return CLIRunner.run_binary(
                settings.TOOL_BINARIES[tool_name], args_list, timeout
            )

        return CLIRunner.run_cli_module(module_path, args_list, timeout)

    @staticmethod
    async def run_cli_module_async(
        module_path: str, args: List[str], timeout: int = None
    ) -> Tuple[int, str, str, List[str]]:
        """Non-blocking wrapper for run_cli_module using ThreadPoolExecutor."""
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            get_cli_executor(),
            lambda: CLIRunner.run_cli_module(module_path, args, timeout)
        )

    @staticmethod
    async def run_cli_with_json_async(
        module_path: str, json_args: str, timeout: int = None
    ) -> Tuple[int, str, str, List[str]]:
        """Non-blocking wrapper for run_cli_with_json using ThreadPoolExecutor."""
        # Resolve tool_name first (outside executor to avoid closure issue)
        tool_name = next((name for name, path in settings.TOOLS_CONFIG if path == module_path), None)

        def convert():
            return CLIRunner._dict_to_args(
                tool_name,
                json.loads(json_args)
            )

        loop = asyncio.get_event_loop()
        args_list = await loop.run_in_executor(get_cli_executor(), convert)

        # Dispatch to binary if this is a binary tool
        if tool_name and tool_name in settings.TOOL_BINARIES:
            return await CLIRunner.run_binary_async(
                settings.TOOL_BINARIES[tool_name], args_list, timeout
            )

        return await CLIRunner.run_cli_module_async(module_path, args_list, timeout)

    @staticmethod
    def run_binary(
        bin_path: str, args: List[str], timeout: int = None
    ) -> Tuple[int, str, str, List[str]]:
        """
        Run an external binary (non-Python) with given args.

        Returns:
            Tuple of (exit_code, stdout, stderr, command)
        """
        if timeout is None:
            timeout = settings.CLI_TOOL_TIMEOUT

        cmd = [bin_path] + args

        try:
            env = CLIRunner._build_env()
            env.pop("PYTHONPATH", None)
            print("\n" + "="*60, flush=True)
            print(f"[BINARY] EXECUTING: {' '.join(cmd)}", flush=True)
            print("="*60 + "\n", flush=True)
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=settings.PROJECT_ROOT,
                env=env,
            )
            print(f"[BINARY] EXIT CODE: {result.returncode}", flush=True)
            if result.stdout:
                print(f"[BINARY] STDOUT:\n{result.stdout}", flush=True)
            if result.stderr:
                print(f"[BINARY] STDERR:\n{result.stderr}", flush=True)
            print("="*60 + "\n", flush=True)
            return result.returncode, result.stdout, result.stderr, cmd
        except subprocess.TimeoutExpired:
            print(f"[BINARY] TIMEOUT: Command timed out after {timeout} seconds", flush=True)
            return -1, "", f"Command timed out after {timeout} seconds", cmd
        except Exception as e:
            print(f"[BINARY] EXCEPTION: {str(e)}", flush=True)
            return -1, "", str(e), cmd

    @staticmethod
    async def run_binary_async(
        bin_path: str, args: List[str], timeout: int = None
    ) -> Tuple[int, str, str, List[str]]:
        """Non-blocking wrapper for run_binary using ThreadPoolExecutor."""
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            get_cli_executor(),
            lambda: CLIRunner.run_binary(bin_path, args, timeout)
        )

    @staticmethod
    def _dict_to_args(tool_name: str, args_dict: Dict[str, Any]) -> List[str]:
        """Convert a dict of arguments to command-line format."""
        args = []
        for key, value in args_dict.items():
            # Determine if the argument is defined as positional or has a short flag
            is_positional = False
            short_flag = None
            if tool_name and tool_name in ARGUMENT_SCHEMAS:
                for arg in ARGUMENT_SCHEMAS[tool_name].get("arguments", []):
                    if arg.get("name") == key:
                        if arg.get("positional"):
                            is_positional = True
                        if "short" in arg:
                            short_flag = arg["short"]
                        break

            if is_positional:
                if isinstance(value, list):
                    for item in value:
                        args.append(str(item))
                else:
                    args.append(str(value))
                continue

            cli_key = short_flag if short_flag else f"--{key.replace('_', '-')}"
            if value is not None:
                if isinstance(value, bool):
                    if value:
                        args.append(cli_key)
                elif isinstance(value, list):
                    if len(value) == 1:
                        args.extend([cli_key, str(value[0])])
                    else:
                        for item in value:
                            args.extend([cli_key, str(item)])
                else:
                    args.extend([cli_key, str(value)])

        return args


# Predefined argument schemas for each CLI tool
# This maps UI forms to command-line arguments

ARGUMENT_SCHEMAS = {
    "bam-basic-stat": {
        "arguments": [
            {"name": "bams", "type": "file", "required": True, "multiple": True, "positional": True},
            {"name": "channel_tag", "type": "select", "default": "ch", "options": ["ch", "zm"]},
            {"name": "min_rq", "type": "number", "required": False, "placeholder": "0.8"},
        ],
    },

    "reads-quality-stats-v3": {
        "arguments": [
            {"name": "bams", "type": "file", "required": True, "multiple": True, "placeholder": "支持通配符 *, 例如：/path/to/*.bam"},
            {"name": "refs", "type": "file", "required": False, "multiple": True, "placeholder": "如果不提供，则不输出对齐相关的指标"},
            {"name": "short_aln", "type": "number", "required": False, "default": 0, "placeholder": "0-200, 用于查询或目标长度"},
            {"name": "force", "type": "boolean", "required": False, "default": False, "help": "重新生成指标文件如果存在"},
        ],
    },
    "sequencing-report-v2": {
        "arguments": [
            {"name": "bams", "type": "file", "required": True, "multiple": True, "placeholder": "支持通配符 *, 例如：/path/to/*.bam"},
            {"name": "refs", "type": "file", "required": False, "multiple": True, "placeholder": "如果不提供，则不输出对齐相关的指标"},
            {"name": "fact_csvs", "type": "file", "required": False, "multiple": True, "placeholder": "合并模式：已有的 fact metric CSV 文件路径"},
            {"name": "basic_csvs", "type": "file", "required": False, "multiple": True, "placeholder": "合并模式：已有的 basic CSV 文件路径"},
            {"name": "short_aln", "type": "number", "required": False, "default": 0, "placeholder": "0-200, 用于查询或目标长度"},
            {"name": "disable_basic_stat", "type": "boolean", "required": False, "default": False, "help": "跳过基础统计计算"},
            {"name": "disable_align_stat", "type": "boolean", "required": False, "default": False, "help": "跳过对齐统计计算"},
            {"name": "np_range", "type": "text", "required": False, "placeholder": "np 范围过滤, 例如: 1:3,5,7:9"},
            {"name": "rq_range", "type": "text", "required": False, "placeholder": "rq 范围过滤, 例如: 0.9:1.1"},
            {"name": "force", "type": "boolean", "required": False, "default": False, "help": "重新生成指标文件如果存在"},
        ],
    },
    "low-q-analysis": {
        "arguments": [
            {"name": "sbr_bam", "type": "file", "required": True},
            {"name": "smc_bam", "type": "file", "required": True},
        ],
    },
    "homo-and-str-ratio": {
        "arguments": [
            {"name": "files", "type": "file", "required": True, "multiple": True, "positional": True, "placeholder": "支持多文件: .bam, .fq, .fastq, .fa, .fasta"},
            {"name": "rq_thr", "type": "number", "required": False, "default": 0.95, "placeholder": "0.0-1.0, 默认 0.95"},
        ],
    },
    "asrtc": {
        "arguments": [
            {"name": "ref_fa", "type": "file", "required": True, "placeholder": "参考序列 FASTA 文件"},
            {"name": "sbr", "type": "file", "required": True, "placeholder": "subreads BAM 文件", "short": "-q"},
            {"name": "smc", "type": "file", "required": True, "placeholder": "SMC BAM/FASTA 文件", "short": "-t"},
            {"name": "prefix", "type": "text", "required": True, "placeholder": "输出前缀 (输出文件: {prefix}.asrtc.txt)", "short": "-p"},
            {"name": "rq_range", "type": "text", "required": False, "placeholder": "rq 范围过滤, 例如: 0.8:1.0"},
            {"name": "np_range", "type": "text", "required": False, "placeholder": "np 范围过滤, 例如: 1:3,5,7:9"},
        ],
    },
}
