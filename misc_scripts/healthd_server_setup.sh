#!/bin/bash

# healthd_service_setup.sh
# 功能：创建并配置 healthd docker 服务的 systemd 自启动

set -e  # 遇到错误时退出脚本


echo "=== HealthD Docker 服务部署脚本 ==="

# 检查是否以 root 或 sudo 权限运行
if [ "$EUID" -ne 0 ]; then
    echo "请使用 sudo 运行此脚本"
    echo "用法: sudo $0"
    exit 1
fi


# 定义变量
SERVICE_NAME="healthd"
SERVICE_FILE="/etc/systemd/system/${SERVICE_NAME}.service"
IMAGE_TAG="192.168.3.38:5000/algo/healthd:1.1.0"
USER_NAME="user"

echo "1. 检查 Docker 服务状态..."
if ! systemctl is-active --quiet docker; then
    echo "错误: Docker 服务未运行，请先启动 Docker"
    exit 1
fi

echo "2. 创建 systemd 服务文件..."
cat > "${SERVICE_FILE}" << EOF
[Unit]
Description=HealthD Docker Service
After=network.target docker.service
Wants=network.target docker.service

[Service]
Type=simple
User=${USER_NAME}
Restart=always
RestartSec=30
ExecStartPre=-/usr/bin/docker rm -f ${SERVICE_NAME}
ExecStart=/usr/bin/docker run \\
    --name ${SERVICE_NAME} \\
    --gpus all \\
    --net=host \\
    -v /proc:/host/proc:ro \\
    -v /sys:/host/sys:ro \\
    -v /:/rootfs:ro \\
    --cap-add=SYS_ADMIN \\
    --privileged \\
    -v /run/udev:/run/udev:ro \\
    -v /var/lib/smartmontools:/var/lib/smartmontools \\
    -v /usr/sbin/smartctl:/usr/sbin/smartctl:ro \\
    -v /home/${USER_NAME}/healthd-log-dir:/healthd-log-dir \\
    -v /home/${USER_NAME}/prometheus-data:/prometheus-data-dir \\
    ${IMAGE_TAG} \\
    start_all

ExecStop=/usr/bin/docker stop -t 10 ${SERVICE_NAME}
ExecStopPost=-/usr/bin/docker rm -f ${SERVICE_NAME}

[Install]
WantedBy=multi-user.target
EOF

echo "服务文件已创建: ${SERVICE_FILE}"


echo "3. 重新加载 systemd 配置..."
systemctl daemon-reload
echo "systemd 配置已重新加载"

echo "4. 启用开机自启动..."
systemctl enable "${SERVICE_NAME}.service"
echo "开机自启动已启用"

echo "5. 启动 ${SERVICE_NAME} 服务..."
if systemctl start "${SERVICE_NAME}.service"; then
    echo "服务启动成功"
else
    echo "警告: 服务启动可能有问题"
fi

echo "6. 检查服务状态..."
sleep 60  # 等待服务完全启动
systemctl status "${SERVICE_NAME}.service" --no-pager

echo ""
echo "=== 部署完成 ==="
echo "服务名称: ${SERVICE_NAME}"
echo "服务文件: ${SERVICE_FILE}"
echo "Docker 镜像: ${IMAGE_TAG}"
echo "运行用户: ${USER_NAME}"
echo ""
echo "常用命令:"
echo "  sudo systemctl status ${SERVICE_NAME}.service    # 查看状态"
echo "  sudo systemctl stop ${SERVICE_NAME}.service      # 停止服务"
echo "  sudo systemctl start ${SERVICE_NAME}.service     # 启动服务"
echo "  sudo systemctl restart ${SERVICE_NAME}.service   # 重启服务"
echo "  sudo journalctl -u ${SERVICE_NAME}.service -f    # 查看日志"
echo "  sudo systemctl disable ${SERVICE_NAME}.service   # 禁用自启动"

# 检查服务是否正常运行
echo ""
echo "检查 Docker 容器状态..."
if docker ps | grep -q "${SERVICE_NAME}"; then
    echo "✓ Docker 容器正在运行"
else
    echo "✗ Docker 容器未运行，请检查日志"
    echo "查看日志命令: sudo journalctl -u ${SERVICE_NAME}.service -n 50"
fi

