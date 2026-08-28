import argparse
import os
import subprocess
import sys


def create_dir_with_777(path: str):
    # 1. 检查目录是否存在
    if not os.path.exists(path):
        # 2. 创建目录 (exist_ok=True 防止并发创建时报错)
        # mode=0o777 在某些系统上受 umask 影响，可能无法直接达到 777
        os.makedirs(path, mode=0o777, exist_ok=True)
        print(f"目录已创建: {path}")
    else:
        print(f"目录已存在: {path}")

    # 3. 显式修改权限为 777
    # 0o777 是八进制表示法
    os.chmod(path, 0o777)
    print(f"权限已设置为 777")


def scp_directory_with_key(
    remote_user,
    remote_host,
    remote_path,
    local_target_path,
    ssh_key_path=None,
    port=22
):
    """
    使用 scp 复制远程目录到本地（使用 SSH 密钥认证）
    :param remote_user: 远程用户名
    :param remote_host: 远程主机 IP 或域名
    :param remote_path: 远程目录路径
    :param local_target_path: 本地目标目录
    :param ssh_key_path: 可选，SSH 私钥路径（如 ~/.ssh/id_rsa）
    :param port: SSH 端口，默认 22
    """
    scp_cmd = ["scp", "-r", "-P", str(port)]

    if ssh_key_path:
        scp_cmd += ["-i", ssh_key_path]

    # 禁用主机密钥检查（自动接受）
    scp_cmd += [
        "-o", "StrictHostKeyChecking=no",
        "-o", "UserKnownHostsFile=/dev/null"
    ]

    create_dir_with_777(local_target_path)

    remote_spec = f"{remote_user}@{remote_host}:{remote_path}"
    scp_cmd += [remote_spec, local_target_path]

    print(f"执行命令: {' '.join(scp_cmd)}")

    try:
        result = subprocess.run(
            scp_cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        print("✅ 目录复制成功！")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ scp 失败:\n{e.stderr}", file=sys.stderr)
        return False


def mkdir_p(path: str):
    """创建目录（含父目录）"""
    os.makedirs(path, exist_ok=True)


def fetch_by_pattern(
    remote_user,
    remote_host,
    base_remote,
    local_target,
    ssh_key_path,
    port,
    pattern,
    full_path=False
):
    """
    根据指定模式拉取远程数据

    :param remote_user: SSH 用户名
    :param remote_host: 远程主机 IP
    :param base_remote: 远程基准目录（如 /data1/EurusResV3/Run0002）
    :param local_target: 本地目标目录
    :param ssh_key_path: SSH 私钥路径
    :param port: SSH 端口
    :param pattern: all=全目录, bam=*.bam, called=*_called.bam, adapter=*_adapter.bam
    :param full_path: True 时 base_remote 为完整远程路径（目录或单个文件），不做前缀拼接且不加尾部斜杠
    """
    # 创建本地目标目录
    mkdir_p(local_target)
    create_dir_with_777(local_target)

    if pattern == "all":
        # === scp -r 拉取整个目录 ===
        base_remote = base_remote.rstrip("/")
        cmd = ["scp", "-r"]

        if ssh_key_path:
            cmd += ["-i", ssh_key_path]
        cmd += ["-P", str(port)]

        # 禁用主机密钥检查
        cmd += [
            "-o", "StrictHostKeyChecking=no",
            "-o", "UserKnownHostsFile=/dev/null"
        ]

        remote_spec = f"{remote_user}@{remote_host}:{base_remote}"
        # full_path 模式下 base_remote 可能是单个文件（不能加尾部斜杠）；
        # 否则视为目录，加尾部斜杠以便 scp 复制到本地目标目录下
        if not full_path:
            remote_spec += "/"
        cmd += [remote_spec, local_target]

        print(f"\n{'='*60}")
        print(f"📥 拉取模式: 全目录 (all)")
        print(f"   远程: {remote_spec}")
        print(f"   本地: {local_target}")
        print(f"   命令: {' '.join(cmd)}")
        print(f"{'='*60}\n")

        try:
            subprocess.run(
                cmd, check=True, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            print("✅ 全目录拉取成功！")
            return True
        except subprocess.CalledProcessError as e:
            print(f"❌ scp 失败:\n{e.stderr}", file=sys.stderr)
            return False

    else:
        # === rsync --include/--exclude 过滤拉取 ===
        # 确定文件名模式
        pattern_map = {
            "bam": "*.bam",
            "called": "*_called.bam",
            "adapter": "*_adapter.bam",
        }
        file_pattern = pattern_map.get(pattern)
        if not file_pattern:
            print(f"❌ 未知的文件模式: {pattern}", file=sys.stderr)
            return False

        base_remote_clean = base_remote.rstrip("/")
        cmd = [
            "rsync", "-avzP",
            f"-i{ssh_key_path}",
            f"-p{port}",
            "--include='*/'",           # 递归进入所有子目录
            f"--include='{file_pattern}'",    # 匹配目标文件
            "--exclude='*'",            # 排除其余
            f"{remote_user}@{remote_host}:{base_remote_clean}/",
            local_target,
            "-o", "StrictHostKeyChecking=no",
            "-o", "UserKnownHostsFile=/dev/null"
        ]

        desc = {"bam": "*.bam", "called": "*_called.bam", "adapter": "*_adapter.bam"}[pattern]
        print(f"\n{'='*60}")
        print(f"📥 拉取模式: {desc}")
        print(f"   远程: {base_remote_clean}/")
        print(f"   本地: {local_target}")
        print(f"   命令: {' '.join(cmd)}")
        print(f"{'='*60}\n")

        try:
            subprocess.run(
                cmd, check=True, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            print("✅ 文件拉取成功！")
            return True
        except subprocess.CalledProcessError as e:
            print(f"❌ rsync 失败:\n{e.stderr}", file=sys.stderr)
            return False


def main():
    """CLI 入口：解析参数并执行数据拉取"""
    parser = argparse.ArgumentParser(
        description="从远端服务器拉取 EurusResV3 数据（支持 SCP/rsync）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 拉取 IP .37 上 Run0002 的全部数据
  python fetch_data_using_scp_cli.py -i 37 -d Run0002

  # 只拉 *.bam 文件
  python fetch_data_using_scp_cli.py -i 37 -d Run0002 -p bam

  # 只拉 *_called.bam 文件
  python fetch_data_using_scp_cli.py -i 37 -d Run0002 -p called

  # 自定义目标目录和 SSH key
  python fetch_data_using_scp_cli.py -i 37 -d Run0002 -t /data1/my_dir -k ~/.ssh/my_key

  # 拉取远程任意完整路径（目录或单个文件，不拼接默认前缀）
  python fetch_data_using_scp_cli.py -i 37 -d /data2/ReProcess/Run0002_called.bam -p all --full-path
"""
    )

    parser.add_argument("--ip", "-i", type=str, required=True,
                        help="IP 地址最后一位，如 37（完整 IP: 192.168.3.<IP>）")
    parser.add_argument("--dir", "-d", type=str, required=True,
                        help="远程目录名，如 Run0002 或 20260722_250302Y0006_Run0002")
    parser.add_argument("--full-path", "--fp", dest="full_path", action="store_true",
                        help="将 --dir/-d 视为远程完整路径（支持任意目录或文件），不拼接 /data1/EurusResV3/ 前缀")
    parser.add_argument("--pattern", "-p", type=str, default="all",
                        choices=["all", "bam", "called", "adapter"],
                        help="拉取模式：all=全目录, bam=*.bam, called=*_called.bam, adapter=*_adapter.bam (默认: all)")
    parser.add_argument("--target", "-t", type=str, default="./data",
                        help="本地目标目录（默认: ./data）")
    parser.add_argument("--user", "-u", type=str, default="user",
                        help="SSH 用户名（默认: user）")
    parser.add_argument("--key", "-k", type=str, default="~/.ssh/id_rsa",
                        help="SSH 私钥路径（默认: ~/.ssh/id_rsa）")
    parser.add_argument("--port", "-P", type=int, default=22,
                        help="SSH 端口（默认: 22）")

    args = parser.parse_args()

    # 拼接远程主机 IP
    remote_host = f"192.168.3.{args.ip}"

    # 构建基准远程路径（--full-path 时视为完整路径，否则拼接默认前缀）
    if args.full_path:
        base_remote = args.dir
    else:
        base_remote = f"/data1/EurusResV3/{args.dir}"

    # 展开 SSH key 路径（处理 ~ 符号）
    ssh_key_path = os.path.expanduser(args.key)

    print(f"\n{'='*60}")
    print(f"🖥️  远程主机: {args.user}@{remote_host}")
    print(f"   远程目录: {base_remote}/")
    print(f"   SSH key:  {ssh_key_path}")
    print(f"   本地目录: {os.path.abspath(args.target)}")
    print(f"{'='*60}\n")

    success = fetch_by_pattern(
        remote_user=args.user,
        remote_host=remote_host,
        base_remote=base_remote,
        local_target=os.path.abspath(args.target),
        ssh_key_path=ssh_key_path,
        port=args.port,
        pattern=args.pattern,
        full_path=args.full_path,
    )

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
