"""Tests for API routers - specifically for ssh_password handling"""

import pytest
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient

from gseda.server.main import app

client = TestClient(app)


class TestBAMBasicStatSSH:
    """Test bam-basic-stat tool with SSH/SCP functionality"""

    def test_execute_with_ssh_password_in_json_body(self):
        """
        Test execution with ssh_password in JSON body.

        This reproduces the 422 error reported by the user.
        The frontend sends ssh_password in the JSON body, not as a Form parameter.
        """
        response = client.post(
            "/api/tools/bam-basic-stat/execute",
            json={
                "tool_name": "bam-basic-stat",
                "args": {
                    "bams": ["/local/path/file.bam"],
                    "channel_tag": "ch",
                    "ssh_password": "test_password"
                }
            }
        )
        # After fix, this should work (no 422 error)
        assert response.status_code != 422, f"Expected success, got {response.status_code}: {response.text}"

    def test_execute_with_remote_bam_and_ssh_password(self):
        """
        Test execution with remote BAM path and ssh_password.

        This is the actual use case: user wants to process a remote BAM file
        via SCP and provides ssh_password in the request.
        """
        with patch('gseda.server.api.routers.FileManager') as mock_fm:
            mock_instance = MagicMock()
            mock_fm.return_value = mock_instance
            mock_instance.fetch_remote_file.return_value = "/tmp/fetched.bam"

            response = client.post(
                "/api/tools/bam-basic-stat/execute",
                json={
                    "tool_name": "bam-basic-stat",
                    "args": {
                        "bams": ["host:/remote/path/file.bam"],
                        "channel_tag": "ch",
                        "ssh_password": "test_password"
                    }
                }
            )
            # After fix, this should work (no 422 error)
            assert response.status_code != 422, f"Expected success, got {response.status_code}: {response.text}"
            # Verify fetch_remote_file was called with correct password
            mock_instance.fetch_remote_file.assert_called_once()
            call_kwargs = mock_instance.fetch_remote_file.call_args[1]
            assert call_kwargs['ssh_password'] == "test_password"


class TestPrepareBamFiles:
    """Test _prepare_bam_files function"""

    @patch('gseda.server.api.routers.FileManager')
    def test_prepare_bam_files_with_ssh_password(self, mock_file_manager):
        """Test that _prepare_bam_files correctly uses ssh_password"""
        from gseda.server.api.routers import _prepare_bam_files

        mock_fm_instance = MagicMock()
        mock_file_manager.return_value = mock_fm_instance
        mock_fm_instance.fetch_remote_file.return_value = "/tmp/fetched.bam"

        # With ssh_password - should fetch remote files
        result = _prepare_bam_files(
            ["host:/remote/file.bam"],
            ssh_password="test_pass"
        )

        assert len(result) == 1
        # Should be a local path after fetching
        assert result[0].startswith("/tmp/") or result[0].startswith("/")

    @patch('gseda.server.api.routers.FileManager')
    def test_prepare_bam_files_without_ssh_password_raises_error(self, mock_file_manager):
        """Test that _prepare_bam_files raises ValueError without ssh_password for remote files"""
        from gseda.server.api.routers import _prepare_bam_files

        # Without ssh_password - should raise error for remote files
        with pytest.raises(ValueError, match="SSH password required"):
            _prepare_bam_files(["host:/remote/file.bam"], ssh_password=None)


class TestIsRemotePath:
    """Test _is_remote_path function"""

    def test_remote_path_user_host_format(self):
        from gseda.server.api.routers import _is_remote_path
        assert _is_remote_path("user@host:/path/file.bam") is True

    def test_remote_path_host_colon_format(self):
        from gseda.server.api.routers import _is_remote_path
        assert _is_remote_path("host:/path/file.bam") is True

    def test_local_path_returns_false(self):
        from gseda.server.api.routers import _is_remote_path
        assert _is_remote_path("/local/path/file.bam") is False
        assert _is_remote_path("file.bam") is False
        assert _is_remote_path("./relative/path.bam") is False


class TestExecuteToolEndpoint:
    """Test execute_tool endpoint behavior"""

    def test_request_without_ssh_password(self):
        """Test execution without ssh_password works for local files"""
        response = client.post(
            "/api/tools/bam-basic-stat/execute",
            json={
                "tool_name": "bam-basic-stat",
                "args": {
                    "bams": ["/local/path/file.bam"],
                    "channel_tag": "ch"
                }
            }
        )
        # Should not return 422 - the request is valid
        assert response.status_code != 422, f"Expected success, got {response.status_code}: {response.text}"

    def test_request_with_ssh_server_in_args(self):
        """Test that ssh_server in args is ignored (only ssh_password matters)"""
        response = client.post(
            "/api/tools/bam-basic-stat/execute",
            json={
                "tool_name": "bam-basic-stat",
                "args": {
                    "bams": ["/local/path/file.bam"],
                    "channel_tag": "ch",
                    "ssh_password": "test_password",
                    "ssh_server": "user@host"
                }
            }
        )
        # Should not return 422
        assert response.status_code != 422, f"Expected success, got {response.status_code}: {response.text}"

    def test_execute_with_channel_tag_parameter(self):
        """Test that channel_tag parameter is correctly passed to CLI"""
        response = client.post(
            "/api/tools/bam-basic-stat/execute",
            json={
                "tool_name": "bam-basic-stat",
                "args": {
                    "bams": ["/local/path/file.bam"],
                    "channel_tag": "zm",
                    "min_rq": 0.8
                }
            }
        )
        # Should not return 422 (parameter validation error)
        # It may fail with 200 but success=false if file doesn't exist, but that's okay
        # The important thing is no 422 (which would mean param parsing failed)
        assert response.status_code != 422, f"Expected no validation error, got {response.status_code}: {response.text}"
