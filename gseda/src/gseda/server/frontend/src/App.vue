<template>
  <div id="app">
    <el-container class="app-container">
      <!-- Sidebar -->
      <el-aside width="280px" class="sidebar">
        <div class="sidebar-header">
          <h2>GSEDA 工具集</h2>
        </div>
        <div class="search-box">
          <el-input
            v-model="searchQuery"
            placeholder="搜索工具..."
            clearable
          >
            <template #prefix>
              <el-icon><Search /></el-icon>
            </template>
          </el-input>
        </div>
        <div class="tools-list">
          <div
            v-for="tool in filteredTools"
            :key="tool.name"
            :class="['tool-item', { active: selectedTool?.name === tool.name }]"
            @click="selectTool(tool)"
          >
            <span>{{ tool.name }}</span>
          </div>
        </div>
      </el-aside>

      <!-- Main Content -->
      <el-main class="main-content">
        <!-- No Tool Selected -->
        <div v-if="!selectedTool" class="empty-state">
          <el-empty description="请选择一个工具开始使用" :image-size="120">
            <template #description>
              <p style="color: #909399; font-size: 14px;">
                从左侧列表选择一个工具，即可开始进行数据分析
              </p>
            </template>
          </el-empty>
        </div>

        <!-- Tool Form -->
        <div v-else class="tool-panel">
          <div class="panel-header">
            <h1>{{ selectedTool.name }}</h1>
            <el-tag size="small" type="info">{{ selectedTool.module_path }}</el-tag>
          </div>

          <!-- Tool-specific form -->
          <component
            :is="getFormComponent()"
            v-if="selectedTool"
            :selectedTool="selectedTool"
            @execute="onExecute"
          ></component>

          <!-- Results -->
          <ResultViewer
            v-if="result"
            :result="result"
            @copy="copyResult"
          />
        </div>
      </el-main>
    </el-container>
  </div>
</template>

<script setup lang="ts">
import { ref, computed } from 'vue'
import { useRouter } from 'vue-router'
import * as ElementPlusIconsVue from '@element-plus/icons-vue'
import axios from 'axios'

// Components
import ToolForm from './components/ToolForm.vue'
import ResultViewer from './components/ResultViewer.vue'

// Types
interface Tool {
  name: string
  module_path: string
}

interface Result {
  success: boolean
  tool_name: string
  stdout: string
  stderr: string
  exit_code: number
}

const router = useRouter()
const tools = ref<Tool[]>([])
const searchQuery = ref('')
const selectedTool = ref<Tool | null>(null)
const result = ref<Result | null>(null)

// Load tools
const loadTools = async () => {
  try {
    const response = await axios.get('/api/tools/list')
    if (response.data.success) {
      tools.value = response.data.data.tools
    }
  } catch (error) {
    console.error('Failed to load tools:', error)
  }
}

// Filter tools by search query
const filteredTools = computed(() => {
  if (!searchQuery.value) return tools.value
  return tools.value.filter(t => t.name.includes(searchQuery.value))
})

// Select a tool
const selectTool = (tool: Tool) => {
  selectedTool.value = tool
  result.value = null
}

// Get form component for current tool
const getFormComponent = () => {
  return ToolForm
}

// Handle tool execution
const onExecute = async (data: any) => {
  if (!selectedTool.value) return

  try {
    const response = await axios.post(
      `/api/tools/${selectedTool.value.name}/execute`,
      data
    )

    result.value = {
      success: response.data.success,
      tool_name: selectedTool.value.name,
      stdout: response.data.data?.stdout || '',
      stderr: response.data.data?.stderr || '',
      exit_code: response.data.data?.exit_code || -1
    }

    if (response.data.success) {
      ElMessage.success('工具执行成功')
    } else {
      ElMessage.error('工具执行失败，请查看错误信息')
    }
  } catch (error: any) {
    result.value = {
      success: false,
      tool_name: selectedTool.value.name,
      stdout: '',
      stderr: error.message || '请求失败',
      exit_code: -1
    }
    ElMessage.error('执行出错：' + (error.message || '未知错误'))
  }
}

// Copy result to clipboard
const copyResult = () => {
  if (result.value) {
    const text = `工具：${result.value.tool_name}
状态：${result.value.success ? '成功' : '失败'}
退出码：${result.value.exit_code}

--- STDOUT ---
${result.value.stdout || '(无输出)'}

--- STDERR ---
${result.value.stderr || '(无错误)'}`

    navigator.clipboard.writeText(text).then(() => {
      ElMessage.success('结果已复制到剪贴板')
    })
  }
}

// Initialize
loadTools()
</script>

<style scoped>
.app-container {
  min-height: 100vh;
}

.sidebar {
  background: #fff;
  border-right: 1px solid #e4e7ed;
  overflow-y: auto;
}

.sidebar-header {
  padding: 20px;
  border-bottom: 1px solid #e4e7ed;
}

.sidebar-header h2 {
  margin: 0;
  color: #303133;
  font-size: 18px;
}

.search-box {
  padding: 15px 20px;
}

.tools-list {
  max-height: calc(100vh - 180px);
  overflow-y: auto;
}

.tool-item {
  padding: 12px 20px;
  cursor: pointer;
  transition: all 0.3s;
  border-left: 3px solid transparent;
}

.tool-item:hover {
  background: #f5f7fa;
  border-left-color: #409eff;
}

.tool-item.active {
  background: #ecf5ff;
  border-left-color: #409eff;
  color: #409eff;
}

.main-content {
  padding: 20px;
  overflow-y: auto;
}

.empty-state {
  display: flex;
  justify-content: center;
  align-items: center;
  min-height: calc(100vh - 150px);
}

.tool-panel {
  max-width: 1200px;
}

.panel-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 30px;
  padding-bottom: 20px;
  border-bottom: 1px solid #e4e7ed;
}

.panel-header h1 {
  margin: 0;
  font-size: 24px;
  color: #303133;
}
</style>
