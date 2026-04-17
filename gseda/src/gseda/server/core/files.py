"""File Manager - Handle remote and uploaded files for GSEDA Server"""

import os
import uuid
import tempfile
import logging
from pathlib import Path
from typing import Optional, List
from contextlib import contextmanager

import paramiko

# Configure logging for console output
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class FileManager:
    """Manages file operations for remote fetching and uploads"""

    # Class-level storage for uploaded file paths
    _uploaded_files: List[str] = []

    def __init__(self, temp_dir: Optional[str] = None):
        """
        Initialize FileManager.

        Args:
            temp_dir: Directory for temporary files. If None, uses system temp.
        """
        self.temp_dir = temp_dir or tempfile.gettempdir()
        self._ssh_client: Optional[paramiko.SSHClient] = None

    def _generate_temp_path(self, original_name: str) -> str:
        """Generate a unique temp file path with original extension."""
        ext = Path(original_name).suffix if original_name else ""
        unique_name = f"gseda_{uuid.uuid4().hex}{ext}"
        return str(Path(self.temp_dir) / unique_name)

    def fetch_remote_file(
        self, remote_path: str, ssh_password: str, ssh_user: str = "root"
    ) -> str:
        """
        Fetch a file from a remote server via SCP.

        Args:
            remote_path: Remote file path in format 'user@host:/path/to/file'
            ssh_password: SSH password for authentication
            ssh_user: SSH username (default: "root")

        Returns:
            Local path to the downloaded file

        Raises:
            ValueError: If remote_path format is invalid
            paramiko.AuthenticationException: If authentication fails
            Exception: For other SCP errors
        """
        # Parse remote path
        if "@" in remote_path:
            user, rest = remote_path.split("@", 1)
            host, remote_file = rest.split(":", 1)
        else:
            user = ssh_user
            if ":" not in remote_path:
                raise ValueError(
                    "Invalid remote path format. Expected 'user@host:/path' or 'host:/path'"
                )
            host, remote_file = remote_path.split(":", 1)

        local_path = self._generate_temp_path(Path(remote_file).name)

        # Log SCP connection info (without password)
        logger.info(f"[SCP] Connecting to {host} as {user}")
        logger.info(f"[SCP] Remote file: {remote_file}")
        logger.info(f"[SCP] Local temp path: {local_path}")

        try:
            # Create SSH client
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

            # Connect with password
            logger.info(f"[SCP] Authenticating...")
            ssh.connect(hostname=host, username=user, password=ssh_password)
            logger.info(f"[SCP] Authentication successful")

            # SCP transfer
            logger.info(f"[SCP] Starting file transfer...")
            sftp = ssh.open_sftp()
            sftp.get(remote_file, local_path)
            sftp.close()
            ssh.close()
            logger.info(f"[SCP] File transfer completed successfully")

            # Track for cleanup
            FileManager._uploaded_files.append(local_path)

            return local_path

        except paramiko.AuthenticationException:
            raise ValueError(
                f"Authentication failed for {host}. Please check your credentials."
            )
        except Exception as e:
            raise Exception(f"Failed to fetch remote file: {str(e)}")

    def handle_uploaded_file(self, file_path: str, original_name: str) -> str:
        """
        Handle an uploaded file by moving it to temp directory.

        Args:
            file_path: Path to the uploaded file
            original_name: Original filename from upload

        Returns:
            New path to the file in temp directory
        """
        import shutil

        local_path = self._generate_temp_path(original_name)

        # Move file to temp dir
        shutil.move(file_path, local_path)

        # Track for cleanup
        FileManager._uploaded_files.append(local_path)

        return local_path

    @classmethod
    def cleanup_temp_files(cls) -> int:
        """
        Clean up all tracked temporary files.

        Returns:
            Number of files cleaned up
        """
        count = 0
        for file_path in cls._uploaded_files:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    count += 1
            except Exception:
                # Skip files that can't be deleted
                pass

        cls._uploaded_files.clear()
        return count

    @classmethod
    @contextmanager
    def temporary_file(cls, file_path: str):
        """
        Context manager for handling temporary files.

        Args:
            file_path: Path to file that should be cleaned up after use
        """
        try:
            yield file_path
        finally:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
            except Exception:
                pass

    def cleanup(self):
        """Instance method to cleanup tracked files."""
        return self.cleanup_temp_files()


# Global instance for use in FastAPI context
_file_manager: Optional[FileManager] = None


def get_file_manager() -> FileManager:
    """Get or create global FileManager instance."""
    global _file_manager
    if _file_manager is None:
        _file_manager = FileManager()
    return _file_manager


def reset_file_manager():
    """Reset global file manager (useful for testing)."""
    global _file_manager
    if _file_manager:
        _file_manager.cleanup()
    _file_manager = None
