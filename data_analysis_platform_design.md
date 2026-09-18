# 内部 AI 数据分析平台需求文档

**版本：V1.3**（2026-09-17：Platform Server 已实现并通过单元测试；Docker 构建 / 端到端联调待 Docker 环境）

---

# 0. 当前项目现状与实现状态（截至 2026-09-17）

> 本文档 §1 之后的内容是 **目标设计**（要建设的最终形态）。本节先说明**当前仓库实际处于什么状态**，
> 用来区分"已经有的"和"还没做的"。

## 0.1 仓库现状

- 仓库名为 **`gsda`**（不是本文档早期版本里假设的 `analysis-agent/` 全新仓库）。它目前是一个
  **生物信息学 SMC / barcode 分析仓库**，本文档描述的平台能力已在这个仓库之上**开始实现**。
- Platform Server 已落地在 **`gsda_platform/`**（⚠️ 包名不是文档 §47 写的 `platform/`——
  顶层 `platform/` 会遮蔽 Python 标准库 `platform` 模块，导致 `import fastapi` 失败，故改名为
  `gsda_platform/`，模块结构不变）。
- `Dockerfile`、`docker-compose.yml`、`docker/entrypoint.sh` 已写好（待有 Docker 环境的机器上
  首次 `docker build` 验证）；`harness/` 放的是 **Harness 启动契约与说明**，不是重新实现的 Agent。

### 当前目录结构（实际）

```text
gsda/
├── .claude/
│   ├── skills/            # ✅ 已存在：分析 Skills（fastq2bam、barcode_ref_align、smc_* 等）
│   └── agents/            # Claude Code agents（非平台 Agent）
├── gsda_platform/         # ✅ 已实现：Platform Server（auth / runtime / proxy / admin）
│   ├── main.py            #    FastAPI 入口（ROLE=platform）
│   ├── config.py          #    环境驱动配置
│   ├── db/                #    SQLAlchemy + SQLite
│   ├── auth/  runtime/  proxy/  admin/
│   └── tests/             #    35 个单元/集成测试（含 fake Docker 后端，无需 Docker）
├── harness/               # ✅ 已建立：Harness(dsh) 启动契约 README + web.config.yml
├── docker/entrypoint.sh   # ✅ 已写：ROLE=platform / harness 双角色分发
├── scripts/               # ✅ 已存在：run_smc*.sh 等驱动脚本
├── third_party/           # ✅ git 子模块：gseda / gsetl / asts / mm2
├── requirements.txt       # ✅ 已补充平台依赖（docker、argon2-cffi 等）
├── Dockerfile             # ✅ 已写（单 Image 双 Role；待 build 验证）
├── docker-compose.yml     # ✅ 已写（Platform 容器 + docker.sock 挂载）
└── README.md
```

## 0.2 各组件实现状态总览

| 组件 | 对应章节 | 当前状态 |
|---|---|---|
| 分析能力（Skills / Scripts / CLI） | §29–31 | ✅ gsda 已有（`.claude/skills/`、`scripts/`、`third_party/*`） |
| Harness（Agent Runtime + Web UI） | §28, §32–35 | ✅ 采用现成 `DeepSeek Harness`(`dsh`)，`harness/` 存启动契约；容器内由 `dsh web` 启动 |
| Platform Server（Auth / Runtime / Proxy / Admin） | §7, §12–27, §33–38 | ✅ 已实现于 `gsda_platform/`（⚠️ 非 `platform/`，见 §0.1） |
| Runtime Manager（用户→容器生命周期） | §19–27 | ✅ 已实现（create/reuse/crash-recover/idle-stop，§36 并发保护） |
| 单 Image、双 Role（platform / harness） | §5–6, §42–44 | ✅ Dockerfile + entrypoint.sh + compose 已写（待 Docker 环境验证 build/run） |
| 用户数据隔离（Workspace / Results / dsh-home） | §25–27, §51 | ✅ 按内部 ID 挂载；隔离依赖容器边界（待 Docker 联调验证） |
| 反向代理（HTTP / WebSocket / Streaming） | §33–34 | ✅ 已实现（HTTP+WS+流式，保留 public Host 供 dsh 信任围栏） |
| Admin | §38 | ✅ 已实现（users / runtimes / restart / stop / delete） |

## 0.3 关于 "Harness" 的说明（已核实 dsh v0.1.5-rc.2）

- 本文档把用户侧 Agent 统称为 **Harness**。已确定采用现成的 **`DeepSeek Harness`（`dsh`）**
  作为 Agent Runtime（Agent Loop / LLM / Tool Calling / Skills / 文件操作 / Session / Streaming / Web UI）。
- 启动方式为容器内 `dsh web --host 0.0.0.0 --port $HARNESS_PORT --no-open --trusted-host $PUBLIC_AUTHORITY`
  （详见 `harness/README.md`，其中列出了逐条核实的 dsh 契约）。
- **三处与本文档原假设不同的关键点**（`harness/README.md` 有完整说明）：
  1. **dsh 没有 `/health` 端点**（§35 的 `GET /health` 不成立）——就绪探测改为「TCP 连通 + `/` 返回 200」。
  2. **`--trusted-host` 必传**：dsh 的 `/api` 有「浏览器信任围栏」按请求 Host 校验；反向代理后
     Host 是平台公网域名，必须通过 `--trusted-host` 告知 dsh，且代理不得改写 Host，否则 `/api` 全被拒。
  3. 默认绑定 `127.0.0.1:3080`，容器内必须 `--host 0.0.0.0` 才能被平台访问。
- 用户状态（`profiles/`、`sessions/`、`storages/`、`settings.yaml`）位于 `DSH_HOME`，
  Runtime Manager 将其挂到用户持久目录 `{DATA_DIR}/users/<id>/dsh-home`，跨重启/镜像升级保留。

---

# 1. 项目背景

公司已经积累了一批数据分析能力，包括：

* 数据分析脚本
* Python / Rust 工具
* 生物信息学工具
* Claude Code Skills
* 分析 Workflow
* 领域知识
* 数据处理和结果生成工具

目前这些能力主要通过 Claude Code 使用。

项目目标是将这些能力封装成一个公司内部的 **AI Data Analysis Platform**，让公司员工通过浏览器直接使用这些能力。

用户无需：

* 安装 Claude Code
* 安装 Python / Rust 环境
* 配置分析工具
* 配置 Skills
* 手动执行复杂脚本
* 学习命令行

用户只需要：

```text
注册 → 登录 → 打开 Agent → 描述需求 → 获得分析结果
```

系统采用 **DeepSeek Harness** 作为 Agent Runtime 和主要 Web UI。

---

# 2. 项目目标

平台需要解决三个核心问题：

## 2.1 能力共享

将现有：

```text
Skills
Scripts
Tools
Workflow
```

封装进统一的 Agent Image。

---

## 2.2 多用户隔离

每个用户拥有：

```text
独立 Runtime
独立 Container
独立 Workspace
独立 Session
独立分析结果
```

不同用户之间不能互相访问数据。

---

## 2.3 简单部署和维护

平台本身希望：

> **只维护一个 Git Repository 和一个主要 Docker Image。**

不希望为了 Platform Server 和 Agent Runtime 维护两套完全独立的项目。

---

# 3. 总体设计

系统由两个逻辑层组成：

```text
                    Internal AI Platform
                           │
            ┌──────────────┴──────────────┐
            │                             │
      Platform Server                User Runtime
            │                             │
      用户/Runtime管理              DeepSeek Harness
            │                             │
            │                         Skills
            │                         Scripts
            │                         Tools
            │                         Workflow
            │                             │
            └──────── Docker ─────────────┘
```

---

# 4. 单一 Repository 设计

Platform Server 和 Harness 不拆成两个 Git 项目，统一放进现有仓库 `gsda`。

> **原则（已定）：** **不改动现有目录结构，只允许新增目录。**
>
> - `platform/`、`harness/`、`docker/` 是**新增**目录 —— 平台要新建的。
> - `.claude/skills/` **保持原样，不迁移到顶层 `skills/`**：当前项目仍在用 Claude Code，
>   Skills 需继续留在 `.claude/skills/` 才能被 Claude Code 直接加载。打包进 Image 时
>   从 `.claude/skills/` 读取即可，无需重排。
> - `scripts/`、`third_party/`（gseda / gsetl / asts / mm2 子模块）等现有目录**一律不动**。

仓库：

```text
gsda/
```

目标结构（仅标注**新增**的部分，其余维持现状）：

```text
gsda/
│
├── platform/              # 🆕 新增：Platform Server
│   ├── auth/
│   ├── runtime/
│   ├── proxy/
│   ├── admin/
│   ├── db/
│   └── main.py
│
├── harness/               # 🆕 新增：DeepSeek Harness（Agent Runtime + Web UI）
│   └── config/
│
├── docker/                # 🆕 新增
│   ├── entrypoint.sh
│   └── config/
│
├── docker-compose.yml     # 🆕 新增
│
├── .claude/               # ✅ 保持现状：skills/ 不迁移（Claude Code 继续使用）
├── scripts/               # ✅ 保持现状（run_smc*.sh 等，不改动）
├── third_party/           # ✅ 保持现状（gseda / gsetl / asts / mm2 子模块）
├── requirements.txt       # ✅ 保持现状
├── Dockerfile             # 🔧 由空文件填充为多 Role 构建
└── README.md
```

所有平台代码、Skills、Scripts、Harness 配置仍在同一个 Repository 中，但
**不重排现有目录，只叠加新增目录**。

---

# 5. 单一 Docker Image

默认只维护一个主要 Image：

```text
company/analysis-agent:<version>
```

例如：

```text
company/analysis-agent:1.0.0
company/analysis-agent:1.1.0
company/analysis-agent:1.2.0
```

Image 中包含：

```text
Platform Server
DeepSeek Harness
Skills
Scripts
Workflow
分析工具
运行时依赖
```

---

# 6. 一个 Image、两种 Runtime Role

虽然只维护一个 Image，但不能把所有用户运行在同一个 Container 中。

Image 根据启动角色运行不同程序。

```text
                  company/analysis-agent:1.0.0
                              │
               ┌──────────────┴──────────────┐
               │                             │
               ▼                             ▼
        ROLE=platform                  ROLE=harness
               │                             │
               ▼                             ▼
        Platform Server                DeepSeek Harness
```

---

# 7. Platform Container

系统启动时运行一个 Platform Container：

```text
company/analysis-agent:1.0.0
          │
          ▼
     Platform C
          │
          └── Platform Server
```

Platform Server 负责：

* 用户注册
* 用户登录
* Session
* Runtime Manager
* Docker Container 管理
* Harness Proxy
* 管理员功能

Platform Container 是整个系统的控制平面。

---

# 8. Harness Container

每个用户拥有独立 Harness Container。

例如：

```text
                         Image
                          │
              ┌───────────┴───────────┐
              │                       │
              ▼                       ▼
        Harness Alice             Harness Bob
              │                       │
         DeepSeek Harness        DeepSeek Harness
              │                       │
         Alice Workspace         Bob Workspace
```

所有 Harness Container 使用相同 Image：

```text
company/analysis-agent:1.0.0
```

但是 Container 独立。

---

# 9. 为什么不能把 Platform 和所有 Harness 放进同一个 Container

禁止以下架构：

```text
Platform
   +
Alice Harness
   +
Bob Harness
   +
Charlie Harness
```

全部运行在一个 Container 中。

原因：

* 无法实现有效用户隔离
* Container 崩溃影响所有用户
* Workspace 边界不清晰
* Runtime 生命周期无法独立管理
* 资源限制无法独立配置

正确架构：

```text
                    Platform Container
                           │
                     Docker API
                           │
           ┌───────────────┼────────────────┐
           │               │                │
           ▼               ▼                ▼
       Alice C          Bob C           Charlie C
           │               │                │
       Harness          Harness           Harness
```

---

# 10. Platform Container 的特殊权限

Platform Server 需要创建和管理用户 Container。

因此 Platform Container 需要访问 Docker Engine API。

V1 可以使用：

```text
/var/run/docker.sock
```

挂载到 Platform Container。

例如：

```text
Host
 │
 ├── Docker Engine
 │
 └── /var/run/docker.sock
          │
          ▼
   Platform Container
```

Platform Container 可以：

* 创建 Container
* 启动 Container
* 停止 Container
* 删除 Container
* 查询 Container 状态

---

# 11. Harness Container 的权限限制

Harness Container **禁止访问 Docker Socket**。

禁止：

```text
/var/run/docker.sock
```

禁止：

```text
privileged
```

禁止：

```text
host filesystem
```

禁止：

```text
host network
```

Harness 只能：

* 访问自己的 Workspace
* 执行自己的分析脚本
* 使用自己的工具
* 调用模型
* 生成自己的结果

---

# 12. Authentication

公司没有 SSO。

由于平台是公司内网系统，V1 使用简单的本地账户系统。

用户可以自行注册。

注册不进行：

* 邮箱验证
* 手机验证
* 实名验证
* 管理员审核
* 公司身份验证

唯一业务规则：

> Username 不能重复。

---

# 13. User 数据模型

SQLite：

```sql
CREATE TABLE users (
    id INTEGER PRIMARY KEY,
    username TEXT UNIQUE NOT NULL,
    password_hash TEXT NOT NULL,
    is_admin BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

# 14. 用户注册

接口：

```http
POST /api/auth/register
```

请求：

```json
{
    "username": "alice",
    "password": "xxxxxxxx"
}
```

流程：

```text
Register
   │
   ▼
检查 username
   │
   ├── 已存在 → Error
   │
   └── 不存在
          │
          ▼
     Password Hash
          │
          ▼
      Create User
```

---

# 15. Username 规则

建议：

```text
3~32 characters
```

允许：

```text
a-z
A-Z
0-9
_
-
```

Username 大小写统一处理。

例如：

```text
Alice
alice
ALICE
```

视为同一个用户。

---

# 16. Password

禁止保存明文密码。

使用：

```text
Argon2id
```

或：

```text
bcrypt
```

V1：

```text
password >= 8 characters
```

不做复杂密码策略。

---

# 17. Login

接口：

```http
POST /api/auth/login
```

请求：

```json
{
    "username": "alice",
    "password": "xxxxxxxx"
}
```

成功后创建 Server-side Session。

使用：

```text
HttpOnly Cookie
```

Cookie 至少设置：

```text
HttpOnly
SameSite=Lax
```

---

# 18. Session

数据库：

```sql
CREATE TABLE sessions (
    id TEXT PRIMARY KEY,
    user_id INTEGER NOT NULL,
    created_at TIMESTAMP NOT NULL,
    expires_at TIMESTAMP NOT NULL,
    last_active TIMESTAMP NOT NULL
);
```

关系：

```text
Session
   │
   └── User
```

退出：

```http
POST /api/auth/logout
```

删除 Session。

---

# 19. Runtime Manager

Runtime Manager 是平台的核心模块。

负责：

```text
User
 ↓
Runtime
 ↓
Container
 ↓
Harness
```

---

# 20. Runtime 数据模型

```sql
CREATE TABLE runtimes (
    id INTEGER PRIMARY KEY,
    user_id INTEGER UNIQUE NOT NULL,
    container_id TEXT,
    image TEXT NOT NULL,
    workspace_path TEXT NOT NULL,
    results_path TEXT NOT NULL,
    status TEXT,
    created_at TIMESTAMP,
    last_active TIMESTAMP
);
```

V1：

> 一个用户最多一个长期 Runtime。

未来可以扩展为：

```text
User
 ├── Runtime A
 ├── Runtime B
 └── Runtime C
```

---

# 21. Runtime 创建

用户第一次进入 Agent：

```text
Login
  ↓
GET /api/runtime
  ↓
Runtime Manager
  ↓
Runtime 不存在
  ↓
创建 Workspace
  ↓
docker create
  ↓
docker start
  ↓
Health Check
  ↓
Harness Ready
  ↓
返回 Runtime
```

---

# 22. Runtime 复用

如果 Runtime 已经存在：

```text
User
 ↓
Runtime Manager
 ↓
Container exists?
```

如果：

```text
RUNNING
```

直接复用。

如果：

```text
STOPPED
```

重新启动。

---

# 23. Runtime 状态

支持：

```text
CREATING
RUNNING
IDLE
STOPPED
CRASHED
FAILED
```

---

# 24. Runtime 自动回收

V1：

```text
idle > 30 min
```

自动：

```text
stop container
```

但是：

```text
Workspace
Results
```

不删除。

用户下一次访问：

```text
STOPPED
   ↓
docker start
   ↓
Harness Ready
   ↓
RUNNING
```

---

# 25. Workspace

用户 Workspace 使用用户内部 ID，而不是 username。

例如：

```text
/data/users/1001/workspace
/data/users/1001/results
```

而不是：

```text
/data/users/alice
```

这样未来允许修改 username。

---

# 26. Container Volume

Alice：

```text
/data/users/1001/workspace
            ↓
       /workspace

/data/users/1001/results
            ↓
       /results
```

Bob：

```text
/data/users/1002/workspace
            ↓
       /workspace

/data/users/1002/results
            ↓
       /results
```

---

# 27. Image 与用户数据分离

Image：

```text
/opt/analysis
```

保存：

```text
Skills
Scripts
Tools
Dependencies
Harness
```

用户数据：

```text
/workspace
/results
```

保存：

```text
用户输入
分析结果
中间文件
报告
```

Image 不保存用户数据。

---

# 28. Harness Backend

每一个 Runtime 内运行：

```text
DeepSeek Harness
```

Harness 负责：

* Agent Loop
* LLM 调用
* Tool Calling
* Skills
* 文件操作
* Session
* Streaming
* Web UI

Platform Server 不重新实现这些功能。

---

# 29. Skills

- **仓库内（源）：** `.claude/skills/` —— **保持现状，不迁移**（§4 原则）：
  当前项目仍在用 Claude Code，Skills 需留在 `.claude/skills/` 才能被直接加载。
- **Image 内（打包后）：** 构建时从 `.claude/skills/` 复制进 Image 的 `/opt/analysis/skills/`。

仓库内现有 Skills（示例）：

```text
.claude/skills/
├── fastq2bam/
├── barcode_ref_align/
├── smc_barcode_split/
├── smc_ab_test_with_barcode/
├── asrtc_analysis/
└── ...
```

打包进 Image 后：

```text
/opt/analysis/skills/
└── ...（与 .claude/skills/ 一致）
```

Skill 应描述：

```text
Purpose
Input
Workflow
Tools
Output
Error handling
```

---

# 30. Scripts

- **仓库内（源）：** `scripts/` —— **保持现状，不迁移**（§4 原则）。
- **Image 内（打包后）：** 构建时复制进 Image 的 `/opt/analysis/scripts/`。

仓库内现有脚本（示例）：

```text
scripts/
├── run_smicing.sh
├── run_smc.sh
├── run_smc_bystrand.sh
├── run_bam2fx_batch.sh
└── ...
```

脚本应该尽量：

* 参数化
* 不依赖固定绝对路径
* 不使用共享临时目录
* 输出到指定 Job / Workspace

---

# 31. Job 目录

建议分析任务使用：

```text
/workspace/jobs/<job_id>/
```

例如：

```text
/workspace/jobs/20260917-001/
├── input/
├── intermediate/
├── output/
└── report.html
```

这样可以避免多个分析任务相互覆盖。

---

# 32. Harness Web UI

用户只访问统一入口：

```text
https://analysis.company.internal
```

用户不需要知道：

```text
container ID
container IP
workspace path
Harness port
```

---

# 33. Harness Proxy

请求流程：

```text
Browser
   │
   ▼
Platform Server
   │
   ▼
Authentication
   │
   ▼
user_id
   │
   ▼
Runtime Manager
   │
   ▼
container_id
   │
   ▼
Harness
```

Platform Server 负责反向代理到正确的 Harness。

---

# 34. Proxy 必须支持

由于 Harness Web UI 可能使用实时通信，Proxy 必须支持：

```text
HTTP
HTTPS
WebSocket
Streaming
```

不能只实现简单 HTTP 请求转发。

---

# 35. Runtime 与 Harness Ready

> **⚠️ 实现更正（见 `harness/README.md`）：** 实际选用的 `dsh` **没有 `/health` 端点**。
> 因此「Health Check」不是 `GET /health`，而是 **TCP 连通 + `GET /` 返回 HTTP 200**
> （`gsda_platform/runtime/docker.py` 的 `readiness_probe`）。只有该探测成功，
> `Runtime.status` 才置为 `RUNNING`。

Container 启动之后：

```text
Container RUNNING
```

不代表 Harness 已经可以使用。

必须执行：

```text
Health Check
```

例如：

```http
GET /health
```

只有 Harness Ready 后：

```text
Runtime.status = RUNNING
```

浏览器才可以进入。

---

# 36. 并发创建保护

必须防止同一个用户同时创建多个 Runtime。

例如：

```text
Browser A ─┐
            ├── Alice
Browser B ─┘
```

只能产生：

```text
Alice
  ↓
Runtime A
```

不能产生：

```text
Runtime A
Runtime B
```

数据库：

```text
UNIQUE(user_id)
```

同时结合事务 / Lock。

---

# 37. Container 崩溃恢复

如果：

```text
Harness Container
       ↓
Crash
```

Runtime Manager 检测：

```text
CRASHED
```

用户再次访问：

```text
CRASHED
 ↓
restart
 ↓
health check
 ↓
RUNNING
```

超过重试次数：

```text
FAILED
```

---

# 38. Admin

V1 可以使用非常简单的管理员模型：

```sql
users.is_admin
```

管理员可以：

```text
查看用户
查看 Runtime
查看 Container
停止 Runtime
重启 Runtime
删除 Container
切换 Image
```

普通用户不能操作其他用户 Runtime。

---

# 39. Runtime 资源限制

Runtime 可以配置：

```text
CPU
Memory
GPU
Disk
```

例如：

```yaml
resources:
  cpu: 8
  memory: 32G
  gpu: 0
```

GPU Runtime：

```yaml
resources:
  cpu: 16
  memory: 64G
  gpu: 1
```

V1 可以所有用户使用统一资源配置。

---

# 40. Image Version

Image：

```text
company/analysis-agent:1.0.0
```

Runtime 创建时记录：

```text
runtime.image
```

例如：

```text
Alice → 1.0.0
Bob   → 1.0.0
```

升级：

```text
1.1.0
```

可以让新用户使用：

```text
New User → 1.1.0
```

老用户继续：

```text
Alice → 1.0.0
```

这样 Image 升级不会自动破坏正在运行的用户环境。

---

# 41. Image 发布流程

开发者只维护一个 Repository：

```text
git
 │
 ▼
gsda
 │
 ├── Platform
 ├── Harness
 ├── Skills
 ├── Scripts
 └── Workflow
```

构建：

```bash
docker build \
    -t company/analysis-agent:1.1.0 .
```

发布：

```bash
docker push company/analysis-agent:1.1.0
```

---

# 42. Platform Container 启动

使用：

```text
ROLE=platform
```

例如：

```text
company/analysis-agent:1.1.0
```

启动：

```text
Platform Server
```

---

# 43. Harness Container 启动

Runtime Manager 创建：

```text
ROLE=harness
```

例如：

```text
company/analysis-agent:1.1.0
```

启动：

```text
DeepSeek Harness
```

---

# 44. Entry Point

统一：

```text
docker/entrypoint.sh
```

逻辑：

```text
ROLE=platform
    ↓
启动 Platform Server

ROLE=harness
    ↓
启动 DeepSeek Harness

其他
    ↓
Error
```

---

# 45. 推荐技术栈

V1：

```text
Backend:
FastAPI

Database:
SQLite

ORM:
SQLAlchemy

Password:
Argon2id

Session:
Server-side Cookie Session

Container:
Docker SDK

Reverse Proxy:
FastAPI / Nginx / Caddy

Agent (Harness, §0.3):
DeepSeek Harness

Deployment:
Docker Compose
```

---

# 46. 初期不使用 Kubernetes

V1 不需要 Kubernetes。

结构：

```text
Docker Host
│
├── Platform Container
│
├── Alice Harness Container
├── Bob Harness Container
└── Charlie Harness Container
```

当用户数量和资源调度需求增长后，再迁移 Kubernetes。

Runtime Manager 可以抽象成：

```text
RuntimeManager
       │
       ├── DockerBackend
       │
       └── KubernetesBackend
```

---

# 47. Platform Server 目录结构

```text
platform/
├── main.py
│
├── auth/
│   ├── models.py
│   ├── service.py
│   └── routes.py
│
├── runtime/
│   ├── models.py
│   ├── manager.py
│   ├── docker.py
│   └── routes.py
│
├── proxy/
│   └── routes.py
│
├── admin/
│   └── routes.py
│
└── db/
    └── database.py
```

---

# 48. API

## Authentication

```http
POST /api/auth/register
POST /api/auth/login
POST /api/auth/logout
GET  /api/auth/me
```

## Runtime

```http
GET    /api/runtime
POST   /api/runtime
POST   /api/runtime/start
POST   /api/runtime/stop
DELETE /api/runtime
```

## Admin

```http
GET    /api/admin/users
GET    /api/admin/runtimes
POST   /api/admin/runtime/{id}/restart
POST   /api/admin/runtime/{id}/stop
DELETE /api/admin/runtime/{id}
```

---

# 49. 完整用户流程

## 第一次使用

```text
Browser
   ↓
Register
   ↓
Login
   ↓
Platform Server
   ↓
创建 User
   ↓
创建 Workspace
   ↓
Runtime Manager
   ↓
创建 Harness Container
   ↓
Health Check
   ↓
Harness Ready
   ↓
进入 Agent UI
```

---

## 第二次使用

```text
Browser
   ↓
Login
   ↓
Platform Server
   ↓
查找 Runtime
   ↓
RUNNING?
   │
   ├── YES → 直接进入
   │
   └── NO
        ↓
      Start
        ↓
    Health Check
        ↓
      Harness
```

---

# 50. 多用户流程

Alice：

```text
Alice
 ↓
Runtime A
 ↓
Container A
 ↓
Harness A
 ↓
Workspace A
```

Bob：

```text
Bob
 ↓
Runtime B
 ↓
Container B
 ↓
Harness B
 ↓
Workspace B
```

最终：

```text
                         Platform
                            │
                     Docker Engine
                            │
            ┌───────────────┼───────────────┐
            │               │               │
            ▼               ▼               ▼
        Alice C          Bob C          Charlie C
            │               │               │
         Harness         Harness         Harness
            │               │               │
        Alice WS         Bob WS        Charlie WS
```

---

# 51. 数据隔离要求

必须保证：

```text
Alice Agent
    X
Bob Workspace
```

```text
Bob Agent
    X
Alice Workspace
```

并且：

```text
Harness
    X
Docker Socket
```

```text
Harness
    X
Host Filesystem
```

```text
Harness
    X
Other Containers
```

---

# 52. 日志

Platform Server 至少记录：

```text
User registration
Login
Logout
Runtime creation
Runtime start
Runtime stop
Runtime crash
Runtime deletion
Image version
```

例如：

```text
2026-09-17 22:10 alice login
2026-09-17 22:10 runtime created
2026-09-17 22:10 container started
2026-09-17 22:11 harness ready
```

Harness Agent trajectory 和 Conversation 由 Harness 自身负责。

---

# 53. V1 不实现

为了降低开发复杂度，V1 不实现：

```text
SSO
LDAP
OAuth
邮箱验证
手机验证
管理员人工审核
复杂 RBAC
Kubernetes
分布式调度
GPU Scheduler
计费
复杂 Quota
Shared Workspace
团队 Workspace
```

---

# 54. V2 扩展方向

未来可以增加：

```text
SSO
LDAP
RBAC
Kubernetes
GPU Scheduler
Job Queue
Dataset Management
Shared Workspace
Team Workspace
Artifact Management
Quota
Audit
Image Registry
Image 自动升级
多 Runtime
```

---

# 55. MVP 验收标准

## 用户系统

* [ ] 用户可以自行注册
* [ ] Username 重名检查
* [ ] Password Hash
* [ ] 用户可以登录
* [ ] 用户可以退出
* [ ] Session 可以过期
* [ ] 用户不能访问其他用户 Runtime

## Runtime

* [ ] 第一次使用自动创建 Runtime
* [ ] 自动创建 Workspace
* [ ] 自动创建 Harness Container
* [ ] Container 使用指定 Image
* [ ] Runtime 可以复用
* [ ] STOPPED Runtime 可以重新启动
* [ ] Crash Runtime 可以恢复
* [ ] Runtime 可以自动 Idle Stop

## Harness

* [ ] DeepSeek Harness 正常启动
* [ ] Browser 可以访问 Harness Web UI
* [ ] HTTP 正常
* [ ] WebSocket 正常
* [ ] Streaming 正常
* [ ] Skills 正常加载
* [ ] Scripts 可以正常执行
* [ ] 可以生成分析结果

## 隔离

* [ ] Alice 无法访问 Bob Workspace
* [ ] Bob 无法访问 Alice Workspace
* [ ] Harness 无 Docker Socket
* [ ] Harness 无宿主机敏感目录访问权限
* [ ] 用户 Container 独立
* [ ] 用户 Runtime 独立

## 部署

* [ ] Platform 和 Harness 使用同一个 Git Repository
* [ ] Platform 和 Harness 使用同一个主要 Docker Image
* [ ] 可以通过 `ROLE=platform` 启动 Platform
* [ ] 可以通过 `ROLE=harness` 启动 Harness
* [ ] 只需要维护一套 Skills/Scripts/Dependencies
* [ ] 一个 Image 可以创建多个用户 Runtime

---

# 56. 最终系统模型

整个系统最终可以简化成：

```text
                    Git Repository
                       gsda
                        │
                        ▼
              gsda:1.0  (company/analysis-agent:1.0)
                        │
             ┌──────────┴──────────┐
             │                     │
             ▼                     ▼
       ROLE=platform          ROLE=harness
             │                     │
             ▼                     ▼
     Platform Container      User Containers
             │                     │
       Runtime Manager       DeepSeek Harness
             │                     │
       Docker Engine         Skills / Scripts
             │                     │
       ┌─────┼─────┐          Workspace
       │     │     │
       ▼     ▼     ▼
    Alice   Bob  Charlie
```

平台的职责边界最终确定为：

```text
Platform Server
    = User + Auth + Runtime + Isolation + Proxy

DeepSeek Harness
    = Agent + LLM + Skills + Tools + Session + UI

Analysis Image
    = Skills + Scripts + Dependencies + Harness
```

因此整个系统只有一个核心 Repository 和一个主要 Image，但运行时通过不同 Container 实现控制面与用户 Agent Runtime 的隔离。
