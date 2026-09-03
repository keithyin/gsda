# GSEDA Server - 设计文档

## 需求概述

将 gseda/gseda/src/gseda 中的命令行工具服务化，通过 Web 界面和 API 让用户可以调用这些工具。

### 技术栈要求
- **后端**: Python + FastAPI
- **前端**: Vue3
- **部署**: 服务器端代码写在 gseda/gseda/src/gseda/server/ 目录下

## 实现状态

### 已完成的功能

#### 后端核心
- [x] FastAPI 应用入口 (main.py)
- [x] 配置管理 (config.py)
- [x] CLI 执行器 (core/runners.py)
- [x] Pydantic 数据验证 (core/schema.py)
- [x] REST API 端点 (api/routers.py)

#### 前端界面
- [x] Vue3 应用架构
- [x] 工具列表侧边栏（支持搜索）
- [x] BAM Basic Stat 表单组件
- [x] MSA View 表单组件
- [x] 结果展示和复制功能
- [x] CDN 模式 HTML 界面（无需构建即可使用）

#### CLI 工具覆盖
所有 22 个 CLI 工具均已通过 API 暴露：
- `bam-basic-stat` - 完整前端表单支持
- `msa-view` - 完整前端表单支持
- `fastx-basic-stat`, `reads-*`, `sequencing-report-*`, `seq-n-stats-*` 等 - API 可用

### 项目结构

```
gseda/src/gseda/server/
├── __init__.py                    # Package marker
├── config.py                      # 配置管理（工具列表定义）
├── main.py                        # FastAPI 应用入口
├── requirements.txt               # Python 依赖
├── start_server.py                # 启动脚本
├── verify.py                      # 验证脚本
├── README.md                      # 完整文档
├── QUICKSTART.md                  # 快速开始指南
├── IMPLEMENTATION_SUMMARY.md      # 实现总结
├── .gitignore                     # Git 忽略规则
├── .env.example                   # 环境变量示例
├── api/                           # API 路由层
│   ├── __init__.py
│   └── routers.py                 # REST API 定义
├── core/                          # 核心功能层
│   ├── __init__.py
│   ├── runners.py                 # CLI 执行器
│   └── schema.py                  # Pydantic 模型
├── templates/                     # HTML 模板
│   └── index.html                 # Web 界面（CDN 模式）
├── static/                        # 静态资源目录
│   └── dist/                      # Vue3 构建输出
└── frontend/                      # Vue3 源代码
    ├── package.json               # npm 依赖
    ├── vite.config.ts             # Vite 配置
    ├── tsconfig.json              # TypeScript 配置
    └── src/
        ├── main.ts                # Vue 入口
        ├── App.vue                # 主组件
        ├── router.ts              # 路由配置
        ├── api.ts                 # API 客户端
        ├── types.ts               # 类型定义
        └── components/            # 组件目录
            ├── ToolForm.vue       # 表单容器
            ├── ResultViewer.vue   # 结果查看器
            └── forms/             # 工具特定表单
                ├── BAMBasicStatForm.vue
                └── MSAViewForm.vue
```

## API 文档

### 端点列表

| 方法 | 路径 | 说明 |
|------|------|------|
| GET | `/api/tools/list` | 获取工具列表 |
| GET | `/api/tools/{name}/info` | 获取工具元数据 |
| POST | `/api/tools/{name}/execute` | 执行工具（JSON body） |
| GET | `/api/tools/{name}/execute` | 执行工具（query params） |

### 示例请求

```bash
# 获取所有可用工具
curl http://localhost:8000/api/tools/list

# 执行 bam-basic-stat
curl -X POST "http://localhost:8000/api/tools/bam-basic-stat/execute" \
  -H "Content-Type: application/json" \
  -d '{"tool_name": "bam-basic-stat", "args": {"bams": ["/path/to/file.bam"], "channel_tag": "ch"}}'

# 通过 URL 参数执行（简写）
curl "http://localhost:8000/api/tools/bam-basic-stat/execute?bams=/path/to/file.bam&channel_tag=ch"
```

## 扩展指南

### 添加新 CLI 工具

1. **编辑 config.py**:
   ```python
   TOOLS_CONFIG.append(("new-tool-name", "gseda.module.path:main_cli"))
   ```

2. **可选 - 创建 API 请求验证**:
   ```python
   # core/schema.py
   class NewToolRequest(ToolRequestBase):
       param1: str
       param2: int
   ```

3. **可选 - 添加前端表单组件**:
   - 在 `frontend/src/components/forms/` 创建新表单
   - 在 `ToolForm.vue` 中添加条件渲染

## 启动方式

### 快速启动（CDN 模式）
```bash
cd gseda/src/gseda/server
uvicorn gseda.server.main:app --reload
```
访问 http://localhost:8000

### 完整启动（含构建）
```bash
cd gseda/src/gseda/server/frontend
npm install && npm run build
cd ..
uvicorn gseda.server.main:app
```

## 下一步规划

- [ ] 完善表单验证反馈
- [ ] 添加文件上传功能
- [ ] 任务队列支持（处理长时间运行任务）
- [ ] API 版本控制
- [ ] 用户认证和权限管理
