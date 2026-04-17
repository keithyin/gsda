"""CLI Module Runners - Execute CLI commands via subprocess"""

import subprocess
import json
import logging
from typing import List, Tuple, Any, Dict
from gseda.server.config import settings

# Configure logging for console output
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class CLIRunner:
    """Execute CLI modules and handle their output"""

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
        args_dict = json.loads(json_args)
        args_list = CLIRunner._dict_to_args(args_dict)
        print(f"[CLI RUNNER] Converted to CLI args: {args_list}")
        return CLIRunner.run_cli_module(module_path, args_list, timeout)

    @staticmethod
    def _dict_to_args(args_dict: Dict[str, Any]) -> List[str]:
        """Convert a dict of arguments to command-line format."""
        args = []
        for key, value in args_dict.items():
            # Special handling for known positional arguments (currently only "bams").
            # Positional arguments are passed without a preceding "--" flag.
            if key == "bams":
                if isinstance(value, list):
                    for item in value:
                        args.append(str(item))
                else:
                    args.append(str(value))
                continue

            # Convert key to CLI format (snake_case to kebab-case)
            # Determine if the argument is defined as positional in ARGUMENT_SCHEMAS.
            is_positional = False
            for schema in ARGUMENT_SCHEMAS.values():
                for arg in schema.get("arguments", []):
                    if arg.get("name") == key and arg.get("positional"):
                        is_positional = True
                        break
                if is_positional:
                    break

            if is_positional:
                # Positional arguments are added without a flag.
                if isinstance(value, list):
                    for item in value:
                        args.append(str(item))
                else:
                    args.append(str(value))
                continue

            # Convert key to CLI format (snake_case to kebab-case)
            cli_key = f"--{key.replace('_', '-') }"
            # Add the key-value pair to args
            if value is not None:
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
    "msa-view": {
        "arguments": [
            {"name": "bam", "type": "file", "required": True},
            {"name": "ref_fasta", "type": "file", "required": True},
            {"name": "ref_name", "type": "text", "required": True},
            {"name": "start", "type": "number", "required": False},
            {"name": "end", "type": "number", "required": False},
            {"name": "o_fasta", "type": "file_output", "required": False},
            {"name": "o_pic", "type": "file_output", "required": False},
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
    "low-q-analysis": {
        "arguments": [
            {"name": "sbr_bam", "type": "file", "required": True},
            {"name": "smc_bam", "type": "file", "required": True},
        ],
    },
}
