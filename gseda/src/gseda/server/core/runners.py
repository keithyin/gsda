"""CLI Module Runners - Execute CLI commands via subprocess"""

import subprocess
import json
from typing import List, Tuple, Any, Dict
from gseda.server.config import settings


class CLIRunner:
    """Execute CLI modules and handle their output"""

    @staticmethod
    def run_cli_module(
        module_path: str, args: List[str], timeout: int = None
    ) -> Tuple[int, str, str]:
        """
        Run a CLI module via subprocess.

        Args:
            module_path: Python module path (e.g., "gseda.bam_ana.bam_basic_stat:main_cli")
            args: Command line arguments to pass
            timeout: Execution timeout in seconds (default from settings)

        Returns:
            Tuple of (exit_code, stdout, stderr)
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
        # Build command to import and execute the target callable.
        # Use a simple script that imports the callable and runs it directly.
        # The previous implementation used double braces which produced a literal
        # `{callable_name}` in the executed code, leading to a NameError.
        cmd = [
            "python",
            "-c",
            f"import sys; from {module_name} import {callable_name}; {callable_name}()",
        ] + args

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=settings.PROJECT_ROOT,
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return -1, "", f"Command timed out after {timeout} seconds"
        except Exception as e:
            return -1, "", str(e)

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
    ) -> Tuple[int, str, str]:
        """
        Run a CLI module with JSON-formatted arguments.

        Args:
            module_path: Python module path
            json_args: JSON string containing argument key-value pairs
            timeout: Execution timeout in seconds

        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        args_dict = json.loads(json_args)
        args_list = CLIRunner._dict_to_args(args_dict)
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
}
