import semver
import subprocess
import logging

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y/%m/%d %H:%M:%S",
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def install_bin(cmd):
    logging.info(f"run command: {cmd}")
    inner_status, inner_output = subprocess.getstatusoutput(cmd)
    assert inner_status == 0, f"run command error. {cmd}, error msg:{inner_output}"


def check_version(bin_name: str, expected_version):
    status, output = subprocess.getstatusoutput(f"{bin_name} -V")
    if status != 0:
        return False
    oup = output.strip()
    version_str = oup.rsplit(" ", maxsplit=1)[1]

    logging.info(f"{bin_name} Version: {version_str}")
    version = semver.Version.parse(version_str)
    if version < expected_version:
        return False

    return True


def check_and_install(bin_name: str, expected_version: semver.Version, install_cmd: str):
    if check_version(bin_name, expected_version):
        return
    install_bin(install_cmd)
    assert check_version(
        bin_name, expected_version), f"check {bin_name} version failed"
