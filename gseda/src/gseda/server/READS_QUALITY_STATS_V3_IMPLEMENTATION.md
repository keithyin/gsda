# Reads Quality Stats V3 实现完成

## 实现概述

基于bam-basic-stat的成功经验，已经完成了reads-quality-stats-v3工具的前后端实现。

## 实现内容

### 1. 后端配置

**文件**: `src/gseda/server/core/runners.py`

添加了ARGUMENT_SCHEMAS定义：

```python
"reads-quality-stats-v3": {
    "arguments": [
        {"name": "bams", "type": "file", "required": True, "multiple": True, "placeholder": "支持通配符 *, 例如: /path/to/*.bam"},
        {"name": "refs", "type": "file", "required": False, "multiple": True, "placeholder": "如果不提供，则不输出对齐相关的指标"},
        {"name": "short_aln", "type": "number", "required": False, "default": 0, "placeholder": "0-200, 用于查询或目标长度"},
        {"name": "force", "type": "boolean", "required": False, "default": False, "help": "重新生成指标文件如果存在"},
    ],
}
```

### 2. 前端表单组件

**文件**: `src/gseda/server/frontend/src/components/forms/ReadsQualityStatsV3Form.vue`

实现了完整的表单UI，包含：
- ✅ 本地文件和远程文件(SCP)切换
- ✅ SSH配置支持
- ✅ BAM文件输入（支持通配符）
- ✅ Reference FASTA文件输入（可选）
- ✅ Short alignment长度设置（可选）
- ✅ 强制重新生成选项
- ✅ 文件上传功能

### 3. 工具表单注册

**文件**: `src/gseda/server/frontend/src/components/ToolForm.vue`

注册了ReadsQualityStatsV3Form组件，使其能在工具列表中被选中。

## 功能特性

### 参数说明

1. **BAM Files** (必需)
   - 支持单个文件路径
   - 支持通配符（如: `/path/to/*.bam`）
   - 支持多文件
   - 支持本地文件和远程文件（SCP）

2. **Reference FASTA** (可选)
   - 用于生成对齐相关的指标
   - 如果不提供，则跳过对齐指标输出

3. **Short Alignment** (可选)
   - 范围: 0-200
   - 默认: 0
   - 用于查询或目标在[30, 200]范围内的reads分析

4. **Force Regenerate** (可选)
   - 勾选后强制重新生成指标文件
   - 默认跳过已存在的指标文件

### SSH远程文件支持

- ✅ 完整的SCP支持
- ✅ SSH服务器配置
- ✅ SSH密码安全传输（在request body中，不在URL中）

## 测试步骤

### 1. 重启服务器

```bash
# 停止当前服务器 (Ctrl+C)
cd /root/projects/gsda/gseda/src
uvicorn gseda.server.main:app --reload --host 0.0.0.0 --port 8000
```

### 2. 测试本地文件

1. 打开浏览器访问 `http://localhost:8000`
2. 在左侧工具列表中找到 "reads-quality-stats-v3"
3. 点击选择该工具
4. 填写以下信息：
   - BAM文件路径: `/path/to/test.bam` 或使用通配符 `/data/*.bam`
   - Reference FASTA: `/path/to/reference.fa`（可选）
   - Short Alignment: `0` 或 `100`
   - 勾选 "强制重新生成指标文件"（可选）
5. 点击 "开始分析"

### 3. 测试远程文件（SCP）

1. 切换到 "远程文件 (SCP)" 模式
2. 填写SSH配置：
   - SSH服务器: `user@host`
   - SSH密码: `your_password`
3. 填写BAM文件路径（远程路径）
4. 填写Reference文件路径（远程路径，可选）
5. 设置Short Alignment长度
6. 点击 "开始分析"

### 4. 测试文件上传

1. 切换到 "本地文件" 模式
2. 点击上传区域选择一个或多个BAM文件
3. 文件上传后会自动添加到输入框
4. 填写其他参数
5. 点击 "开始分析"

## 日志检查

测试时请检查：

1. **浏览器Console日志**：
   - 应该看到表单提交的完整参数
   - 应该看到 `[Axios Request]` 日志

2. **服务器Console日志**：
   - 应该看到 `=== routers.py - execute_tool (POST) START ===`
   - 应该看到完整的args参数
   - 应该看到SSH配置信息（如果使用远程文件）

3. **Network面板**：
   - 请求方法应该是 **POST**
   - URL应该是 `http://localhost:8000/api/tools/reads-quality-stats-v3/execute`
   - URL中不应该包含SSH密码
   - Request Body应该包含JSON格式的参数

## 与bam-basic-stat的区别

| 特性 | bam-basic-stat | reads-quality-stats-v3 |
|------|----------------|------------------------|
| 支持通配符 | ❌ | ✅ |
| 需要Reference | ❌ | ✅（可选）|
| Short Alignment参数 | ❌ | ✅ |
| 强制重新生成选项 | ❌ | ✅ |
| 文件上传 | ✅ | ✅ |
| SCP远程文件 | ✅ | ✅ |

## 安全特性

✅ SSH密码不在URL中传输（修复了之前的问题）
✅ 使用HTTPS/HTTP body传输敏感信息
✅ 符合安全最佳实践

## 构建状态

```
✓ frontend构建成功
✓ 所有组件已注册
✓ 参数schema已定义
✓ 无语法错误
```

## 文件清单

### 修改的文件
1. `src/gseda/server/core/runners.py` - 添加ARGUMENT_SCHEMAS
2. `src/gseda/server/frontend/src/components/ToolForm.vue` - 注册新表单组件

### 新增的文件
1. `src/gseda/server/frontend/src/components/forms/ReadsQualityStatsV3Form.vue` - 新工具表单组件

## 注意事项

1. **通配符支持**: 在本地文件输入框中，可以输入通配符如 `*.bam` 或 `/path/*.bam`
2. **Reference文件**: 可选，如果不提供则不输出对齐相关的指标
3. **Short Alignment**: 默认为0，表示不使用short alignment过滤
4. **指标文件**: 如果指标文件已存在，默认会跳过以节省时间。勾选"强制重新生成"可以强制重新计算

## 后续优化建议

1. 添加对多个Reference文件的支持
2. 优化通配符文件匹配的性能
3. 添加进度显示（因为可能处理大量文件）
4. 支持批量文件选择器