# GSEDA Server - 快速开始指南

## 最小化安装

### 1. 后端（立即启动）

```bash
cd gseda/src/gseda/server
pip install fastapi uvicorn jinja2 pydantic python-multipart
uvicorn gseda.server.main:app --reload
```

访问：http://localhost:8000

Web 界面已通过 CDN 加载，无需构建即可使用。

### 2. API 文档

启动服务器后访问：http://localhost:8000/docs

## 完整安装（含前端）

### 1. 安装后端依赖
```bash
cd gseda/src/gseda/server
pip install -r requirements.txt
```

### 2. 安装前端依赖并构建
```bash
cd gseda/src/gseda/server/frontend
npm install
npm run build
```

### 3. 启动服务
```bash
cd gseda/src/gseda/server
uvicorn gseda.server.main:app --host 0.0.0.0 --port 8000
```

## 测试工具执行

### 使用 curl 测试 bam-basic-stat：
```bash
curl -X POST "http://localhost:8000/api/tools/bam-basic-stat/execute" \
  -H "Content-Type: application/json" \
  -d '{
    "tool_name": "bam-basic-stat",
    "args": {
      "bams": ["/path/to/test.bam"],
      "channel_tag": "ch",
      "min_rq": null
    }
  }'
```

### 使用 Python 测试：
```python
import requests

response = requests.post(
    "http://localhost:8000/api/tools/bam-basic-stat/execute",
    json={
        "tool_name": "bam-basic-stat",
        "args": {
            "bams": ["/path/to/test.bam"],
            "channel_tag": "ch"
        }
    }
)
print(response.json())
```

## 添加新的 CLI 工具

只需三步：

1. **编辑 config.py**，在 TOOLS_CONFIG 中添加：
```python
("my-new-tool", "gseda.module.path:main_cli"),
```

2. **可选：添加 API 请求验证**
```python
# core/schema.py
class MyNewToolRequest(ToolRequestBase):
    param1: str
    param2: int
```

3. **可选：添加前端表单组件**
   - 创建 `frontend/src/components/forms/MyNewToolForm.vue`
   - 在 `ToolForm.vue` 中注册该组件

## 常见问题

### Q: Web 界面显示空白？
A: 这是正常的 CDN 模式，需要联网加载 Vue 和 Element Plus。

### Q: CLI 工具执行超时？
A: 在 config.py 中调整 `CLI_TOOL_TIMEOUT` 变量。

### Q: 如何离线使用？
A: 将静态资源下载到 local，替换 CDN 链接为本地路径。

## 下一步

- 阅读 [`README.md`](./README.md) 获取完整文档
- 查看 [`IMPLEMENTATION_SUMMARY.md`](./IMPLEMENTATION_SUMMARY.md) 了解实现细节
- 访问 API 文档：http://localhost:8000/docs
