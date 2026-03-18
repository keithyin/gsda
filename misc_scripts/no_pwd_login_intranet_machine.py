#!/usr/bin/env python3
import os
import sys
import subprocess
import ipaddress
import pexpect # pip install pexpect
import argparse
import time

def run_cmd(cmd, shell=False):
    try:
        return subprocess.run(
            cmd,
            shell=shell,
            capture_output=True,
            text=True,
            timeout=10
        )
    except subprocess.TimeoutExpired:
        return None

def is_host_alive(host, timeout=1):
    """检查主机是否在线（通过 ping）"""
    result = run_cmd(["ping", "-c", "1", "-W", str(timeout), host])
    return result and result.returncode == 0

def read_public_key(key_path):
    """读取公钥内容"""
    if not os.path.exists(key_path):
        raise FileNotFoundError(f"公钥文件不存在: {key_path}")
    with open(key_path) as f:
        return f.read().strip()

def deploy_key_to_host(host, username, password, pub_key, port=22, dry_run=False):
    """
    将公钥部署到目标主机
    """
    if dry_run:
        print(f"[DRY-RUN] 将部署公钥到 {username}@{host}")
        return True

    # 构造远程命令：创建 .ssh 目录并追加公钥
    remote_cmd = (
        "mkdir -p ~/.ssh && "
        "chmod 700 ~/.ssh && "
        f"echo '{pub_key}' >> ~/.ssh/authorized_keys && "
        "chmod 600 ~/.ssh/authorized_keys && "
        "echo 'OK'"
    )

    ssh_cmd = [
        "ssh",
        "-p", str(port),
        "-o", "StrictHostKeyChecking=no",
        "-o", "UserKnownHostsFile=/dev/null",
        "-o", "ConnectTimeout=10",
        f"{username}@{host}",
        remote_cmd
    ]

    try:
        child = pexpect.spawn(" ".join(ssh_cmd), encoding='utf-8', timeout=30)
        i = child.expect(['password:', 'Password:', pexpect.EOF, pexpect.TIMEOUT])
        if i in (0, 1):
            child.sendline(password)
            child.expect(pexpect.EOF, timeout=20)
            output = child.before
            if "OK" in output:
                return True
            else:
                print(f"⚠️ {host}: 命令执行无报错但未返回 OK，可能失败")
                return False
        elif i == 2:
            # 可能已免密，无需密码
            output = child.before
            if "OK" in output:
                return True
            else:
                print(f"⚠️ {host}: 无密码提示，但响应异常: {output[:100]}")
                return False
        else:
            print(f"❌ {host}: 超时或未知错误")
            return False
    except Exception as e:
        print(f"❌ {host}: 异常 - {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="批量部署 SSH 公钥到局域网主机")
    parser.add_argument("--user", "-u", required=True, help="目标主机的用户名")
    parser.add_argument("--password", "-p", required=True, help="目标主机的密码")
    parser.add_argument("--network", "-n", default="192.168.3.0/24", help="局域网网段 (默认: 192.168.3.0/24)")
    parser.add_argument("--port", type=int, default=22, help="SSH 端口 (默认: 22)")
    parser.add_argument("--key", default="~/.ssh/id_rsa.pub", help="本地公钥路径 (默认: ~/.ssh/id_rsa.pub)")
    parser.add_argument("--dry-run", action="store_true", help="仅打印将要执行的操作，不实际部署")

    args = parser.parse_args()

    pub_key_path = os.path.expanduser(args.key)
    try:
        pub_key = read_public_key(pub_key_path)
    except FileNotFoundError as e:
        print(e)
        sys.exit(1)

    print(f"🌐 扫描网段: {args.network}")
    print(f"👤 目标用户: {args.user}")
    if args.dry_run:
        print("📝 模式: DRY-RUN（仅预览）")
    print("-" * 50)

    network = ipaddress.IPv4Network(args.network, strict=False)
    total = 0
    success = 0

    for ip in network.hosts():
        host = str(ip)
        if not is_host_alive(host):
            continue

        total += 1
        print(f"➡️  正在处理 {host}...", end=" ")
        if deploy_key_to_host(host, args.user, args.password, pub_key, args.port, args.dry_run):
            print("✅ 成功")
            success += 1
        else:
            print("❌ 失败")

        # 避免太快触发 SSH 限制
        time.sleep(0.5)

    print("-" * 50)
    print(f"📊 总结: 共发现 {total} 台在线主机，成功配置 {success} 台。")

if __name__ == "__main__":
    main()