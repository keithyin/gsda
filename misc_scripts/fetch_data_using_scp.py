import subprocess
import sys
import os


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


# 示例用法
if __name__ == "__main__":

    uris = [
        # "user@192.168.3.125:/data1/512k_test/20260120_250302Y0004_Run0002/*.bam",
        "user@192.168.3.125:/data1/EurusResV3/20260120_250302Y0004_Run0001/*.bam",
        # "user@192.168.3.125:/data1/EurusResV3/20260112_250302Y0004_Run0001",
        
    ]

    for uri in uris:
        user, ext = uri.split("@")
        remote_host, remote_path = ext.split(":")

        success = scp_directory_with_key(
            remote_user=user,
            remote_host=remote_host,
            remote_path=remote_path,
            local_target_path="/data1/ccs_data/20260120-huahongDPN",
            ssh_key_path="~/.ssh/id_rsa",  # 可选
            port=22
        )
