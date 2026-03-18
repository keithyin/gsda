import subprocess
import sys

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
        "user@192.168.3.125:/data1/EurusResV3/20251224_250302Y0004_Run0001",
        "user@192.168.3.125:/data1/EurusResV3/20251224_250302Y0004_Run0002",
        "user@192.168.3.125:/data1/EurusResV3/20251226_250302Y0004_Run0002",
        "user@192.168.3.125:/data1/EurusResV3/20251229_250302Y0004_Run0003",
        "user@192.168.3.125:/data1/EurusResV3/20251230_250302Y0004_Run0001",
        "user@192.168.3.125:/data1/EurusResV3/20251231_250302Y0004_Run0001",
        "user@192.168.3.63:/data1/EurusResV3/20260106_250302Y0001_Run0003",
        "user@192.168.3.63:/data1/EurusResV3/20260106_250302Y0001_Run0002",
        "user@192.168.3.63:/data1/EurusResV3/20260105_250302Y0001_Run0005",
        "user@192.168.3.63:/data1/EurusResV3/20260108_250302Y0001_Run0001",
        "user@192.168.3.72:/data1/EurusResV3/20260108_240601Y0088_Run0001",
        
    ]
    
    for uri in uris:
        user, ext = uri.split("@")
        remote_host, remote_path = ext.split(":")
        
        success = scp_directory_with_key(
            remote_user=user,
            remote_host=remote_host,
            remote_path=remote_path,
            local_target_path="/data1/ccs_data/20260109-saisuofei-resplit",
            ssh_key_path="~/.ssh/id_rsa",  # 可选
            port=22
        )