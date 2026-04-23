#!/usr/bin/env python3

import subprocess
import concurrent.futures
import os
import getpass

SUBNET = "192.168.3"
IMAGE = "192.168.3.38:5000/algo/healthd:1.1.0"
REMOTE_SCRIPT = "/tmp/healthd_service_setup.sh"
LOCAL_SCRIPT = "/root/projects/gsda/misc_scripts/healthd_server_setup.sh"
SSH_USER = "user"

PARALLEL = 1  # ⚠️ sudo 场景建议降低并发

SSH_OPTS = [
    "-o", "StrictHostKeyChecking=no",
    "-o", "UserKnownHostsFile=/dev/null",
    "-o", "ConnectTimeout=3",
    "-o", "BatchMode=yes"
]

# 🔑 输入 sudo 密码（只输入一次）
SUDO_PASSWORD = getpass.getpass("请输入 sudo 密码: ")


def run_cmd(cmd, input_data=None, timeout=120):
    try:
        result = subprocess.run(
            cmd,
            input=input_data,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
        return result.returncode == 0, result.stdout.decode(errors="ignore"), result.stderr.decode(errors="ignore")
    except Exception as e:
        return False, "", str(e)


def ping_host(ip):
    return subprocess.run(
        ["ping", "-c", "1", "-W", "1", ip],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0


def scan_hosts():
    print("1. 扫描存活主机...")
    hosts = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=PARALLEL) as executor:
        futures = {
            executor.submit(ping_host, f"{SUBNET}.{i}"): f"{SUBNET}.{i}"
            for i in range(1, 255)
        }

        for future in concurrent.futures.as_completed(futures):
            if future.result():
                hosts.append(futures[future])

    print("发现主机:")
    for h in hosts:
        print(" ", h)
    print()

    return hosts


def ssh_cmd(ip, command, input_data=None):
    # ⚠️ 加 -tt 解决 sudo requiretty 问题
    cmd = ["ssh"] + SSH_OPTS + [f"{SSH_USER}@{ip}", command]
    return run_cmd(cmd, input_data=input_data)


def scp_file(ip):
    cmd = ["scp"] + SSH_OPTS + [LOCAL_SCRIPT,
                                f"{SSH_USER}@{ip}:{REMOTE_SCRIPT}"]
    ok, out, err = run_cmd(cmd)
    return ok


def deploy_host(ip):
    print(f"---- [{ip}] ----")

    # SSH 测试
    print("    :> check ssh ok?")

    ok, _, _ = ssh_cmd(ip, "echo ok")
    if not ok:
        print(f"    :> SSH 不可用，跳过")
        return
    print("    :> check docker ok?")
    # docker 检查
    ok, _, _ = ssh_cmd(ip, "command -v docker >/dev/null 2>&1")
    if not ok:
        print(f"    :> 未安装 docker，跳过")
        return

    # 镜像检查
    print("    :> check image exists?")
    ok, _, _ = ssh_cmd(ip, f"docker image inspect {IMAGE} >/dev/null 2>&1")
    if ok:
        print(f"[{ip}] 已存在镜像，跳过")
        return

    print(f"    :> 开始部署...")

    # scp 脚本
    print("    :> scp file")

    if not scp_file(ip):
        print(f"    :> scp 失败")
        return

    # 🚨 sudo 执行（关键逻辑）
    sudo_cmd = f"sudo -S bash {REMOTE_SCRIPT}"

    print("    :> exec deploy shell script")

    ok, out, err = ssh_cmd(
        ip,
        sudo_cmd,
        input_data=(SUDO_PASSWORD + "\n").encode()
    )

    if ok:
        print(f"    :> 部署完成")
    else:
        if "password" in err.lower():
            print(f"    :>sudo 密码错误")
        else:
            print(f"    :> 部署失败")
            print(f"    :>stderr: {err.strip()}")


def main():
    print("=== HealthD 批量部署 (sudo 密码模式）===")

    if not os.path.exists(LOCAL_SCRIPT):
        print(f"错误: 找不到 {LOCAL_SCRIPT}")
        return

    # hosts = scan_hosts()
    hosts = [ "192.168.3.{}".format(i) for i in range(0, 255)]
    # hosts = ["192.168.3.10"]

    print("2. 开始批量部署...\n")
    for host in hosts:
        deploy_host(host)

    print("\n=== 部署完成 ===")


if __name__ == "__main__":
    main()
