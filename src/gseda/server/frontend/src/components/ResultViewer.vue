<template>
  <div class="result-viewer">
    <!-- Result Header -->
    <div class="result-header">
      <div class="header-left">
        <h3>执行结果</h3>
        <el-tag :type="result.success ? 'success' : 'danger'" size="small" class="status-badge">
          {{ result.success ? '成功' : '失败' }}
        </el-tag>
        <span class="exit-code" style="margin-left: 10px; font-size: 13px;">
          退出码：{{ result.exit_code }}
        </span>
      </div>
      <div class="header-actions">
        <el-button link @click="toggleMinimize">
          <el-icon>
            <component :is="isMinimized ? 'Bottom' : 'Top'" />
          </el-icon>
        </el-button>
        <el-button link @click="$emit('update:visible', false)">
          <el-icon><Close /></el-icon>
        </el-button>
      </div>
    </div>

    <!-- Tabs -->
    <div class="tabs-container" v-if="!isMinimized">
      <el-radio-group v-model="activeTab" class="tabs-group">
        <el-radio-button value="summary" :disabled="!hasAnyOutput">汇总</el-radio-button>
        <el-radio-button value="files" v-if="hasFileOutputs">文件</el-radio-button>
        <el-radio-button value="stdout" :disabled="!hasContent(result.stdout)">STDOUT</el-radio-button>
        <el-radio-button value="stderr" :disabled="!hasContent(result.stderr)">STDERR</el-radio-button>
      </el-radio-group>
    </div>

    <!-- Search Bar -->
    <div class="search-bar" v-if="!isMinimized && activeTab !== 'summary'">
      <el-input
        v-model="searchQuery"
        placeholder="搜索输出内容..."
        clearable
        size="small"
        @input="handleSearch"
      >
        <template #prefix>
          <el-icon><Search /></el-icon>
        </template>
        <template #suffix v-if="searchMatches > 0">
          <span class="search-count">{{ searchMatches }}/{{ searchTotal }}</span>
        </template>
      </el-input>
    </div>

    <!-- Result Content -->
    <div class="result-output" v-if="!isMinimized">
      <!-- Summary Tab -->
      <div v-if="activeTab === 'summary'" class="summary-tab">
        <div class="summary-info">
          <div class="info-row">
            <span class="info-label">工具名称:</span>
            <span class="info-value">{{ result.tool_name }}</span>
          </div>
          <div class="info-row">
            <span class="info-label">执行状态:</span>
            <el-tag :type="result.success ? 'success' : 'danger'" size="small">
              {{ result.success ? '成功' : '失败' }}
            </el-tag>
          </div>
          <div class="info-row">
            <span class="info-label">退出码:</span>
            <span class="info-value">{{ result.exit_code }}</span>
          </div>
        </div>

        <!-- STDOUT Preview -->
        <div v-if="hasContent(result.stdout)" class="preview-section">
          <div class="preview-header">
            <span class="preview-title">STDOUT 预览</span>
            <el-button link @click="switchTab('stdout')">
              查看完整 <el-icon><ArrowRight /></el-icon>
            </el-button>
          </div>
          <pre class="preview-content">{{ formatOutputForPreview(result.stdout) }}</pre>
        </div>

        <!-- STDERR Preview -->
        <div v-if="hasContent(result.stderr)" class="preview-section error-section">
          <div class="preview-header">
            <span class="preview-title">STDERR 预览</span>
            <div class="preview-header-actions">
              <el-button size="small" type="primary" link @click="askAIHelp" :loading="aiLoading">
                <el-icon><QuestionFilled /></el-icon>
                让 AI 分析错误
              </el-button>
              <el-button link @click="switchTab('stderr')">
                查看完整 <el-icon><ArrowRight /></el-icon>
              </el-button>
            </div>
          </div>
          <pre class="preview-content">{{ formatOutputForPreview(result.stderr, result.tool_name) }}</pre>
          <!-- AI Analysis Result -->
          <div v-if="aiAnalysis === null" class="ai-analysis-section ai-analysis-loading">
            <div class="ai-analysis-header">
              <span class="ai-analysis-title">
                <el-icon class="loading-icon"><Loading /></el-icon>
                AI 正在分析错误...
              </span>
              <span class="ai-analysis-status">请稍候，通常需 2-3 分钟</span>
            </div>
          </div>
          <div v-else-if="aiAnalysis" class="ai-analysis-section">
            <div class="ai-analysis-header">
              <span class="ai-analysis-title">
                <el-icon><QuestionFilled /></el-icon>
                AI 分析结果
              </span>
              <el-button size="small" link @click="copyAIAnalysis">
                <el-icon><DocumentCopy /></el-icon>
                复制
              </el-button>
            </div>
            <div v-if="aiAnalysis" class="ai-analysis-content" v-html="renderMarkdown(aiAnalysis)"></div>
          </div>
        </div>

        <div v-if="!hasContent(result.stdout) && !hasContent(result.stderr)" class="empty-output">
          <el-empty description="无输出内容" :image-size="80" />
        </div>
        <div v-if="hasContent(result.stdout) || hasContent(result.stderr)" class="preview-close">
          <el-button link @click="$emit('update:visible', false)">
            <el-icon><Close /></el-icon>
            收起预览
          </el-button>
        </div>
      </div>

      <!-- Files Tab -->
      <div v-if="activeTab === 'files' && result.file_outputs && result.file_outputs.length > 0" class="output-tab">
        <div class="output-controls">
          <div class="control-left">
            <span class="output-title">输出文件 ({{ result.file_outputs.length }})</span>
          </div>
          <div class="control-right">
            <el-button size="small" @click="downloadAllFiles">
              <el-icon><Download /></el-icon>
              全部下载
            </el-button>
          </div>
        </div>
        <div class="output-content">
          <div class="file-list">
            <div v-for="file in result.file_outputs" :key="file.name" class="file-item">
              <span class="file-name">{{ file.name }}</span>
              <el-button size="small" type="primary" link @click="downloadFile(file.download_url, file.filename)">
                <el-icon><Download /></el-icon>
                下载
              </el-button>
            </div>
          </div>
        </div>
      </div>

      <!-- STDOUT Tab -->
      <div v-if="activeTab === 'stdout' && hasContent(result.stdout)" class="output-tab">
        <div class="output-controls">
          <div class="control-left">
            <el-button link @click="toggleExpand('stdout')">
              <el-icon><component :is="getExpandIcon('stdout')" /></el-icon>
              {{ isExpanded.stdout ? '收起' : '展开' }}
            </el-button>
          </div>
          <div class="control-right">
            <el-button size="small" @click="copySection('stdout')">
              <el-icon><DocumentCopy /></el-icon>
              复制
            </el-button>
            <el-button size="small" @click="downloadSection('stdout', 'stdout')">
              <el-icon><Download /></el-icon>
              下载
            </el-button>
          </div>
        </div>
        <div class="output-content">
          <pre ref="stdoutPre" class="output-text terminal-output">
            {{ result.stdout }}
          </pre>
        </div>
      </div>

      <!-- STDERR Tab -->
      <div v-if="activeTab === 'stderr' && hasContent(result.stderr)" class="output-tab">
        <div class="output-controls">
          <div class="control-left">
            <el-button link @click="toggleExpand('stderr')">
              <el-icon><component :is="getExpandIcon('stderr')" /></el-icon>
              {{ isExpanded.stderr ? '收起' : '展开' }}
            </el-button>
          </div>
          <div class="control-right">
            <el-button size="small" type="primary" @click="askAIHelp" :loading="aiLoading">
              <el-icon><QuestionFilled /></el-icon>
              让 AI 分析错误
            </el-button>
            <el-button size="small" @click="copySection('stderr')">
              <el-icon><DocumentCopy /></el-icon>
              复制
            </el-button>
            <el-button size="small" @click="downloadSection('stderr', 'stderr')">
              <el-icon><Download /></el-icon>
              下载
            </el-button>
          </div>
        </div>
        <div class="output-content terminal-output">
          <pre ref="stderrPre" class="output-text">
            {{ stderrWithContext }}
          </pre>
        </div>
        <!-- AI Analysis Result -->
        <div v-if="aiAnalysis === null" class="ai-analysis-section ai-analysis-loading">
          <div class="ai-analysis-header">
            <span class="ai-analysis-title">
              <el-icon class="loading-icon"><Loading /></el-icon>
              AI 正在分析错误...
            </span>
            <span class="ai-analysis-status">请稍候，通常需 2-3 分钟</span>
          </div>
        </div>
        <div v-else-if="aiAnalysis" class="ai-analysis-section">
          <div class="ai-analysis-header">
            <span class="ai-analysis-title">
              <el-icon><QuestionFilled /></el-icon>
              AI 分析结果
            </span>
            <el-button size="small" link @click="copyAIAnalysis">
              <el-icon><DocumentCopy /></el-icon>
              复制
            </el-button>
          </div>
          <div v-if="aiAnalysis" class="ai-analysis-content" v-html="renderMarkdown(aiAnalysis)"></div>
        </div>
      </div>

      <!-- Empty State -->
      <div v-if="(activeTab === 'stdout' && !hasContent(result.stdout)) || (activeTab === 'stderr' && !hasContent(result.stderr))" class="empty-output">
        <el-empty description="无内容可显示" :image-size="80" />
      </div>
    </div>

    <!-- Action Buttons (Always Visible) -->
    <div class="action-buttons" v-if="!isMinimized">
      <el-button size="small" @click="copySection('stdout')" :disabled="!result.stdout">
        <el-icon><DocumentCopy /></el-icon>
        复制 STDOUT
      </el-button>
      <el-button size="small" @click="copySection('stderr')" :disabled="!result.stderr">
        <el-icon><DocumentCopy /></el-icon>
        复制 STDERR
      </el-button>
      <el-button size="small" @click="copyAll">
        <el-icon><DocumentCopy /></el-icon>
        复制全部
      </el-button>
      <el-button size="small" @click="downloadAll">
        <el-icon><Download /></el-icon>
        下载全部
      </el-button>
      <el-button size="small" @click="$emit('update:visible', false)">
        <el-icon><Close /></el-icon>
        关闭
      </el-button>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, computed, watch, nextTick } from 'vue'
import { ElMessage } from 'element-plus'
import { marked } from 'marked'
import {
  DocumentCopy,
  Download,
  Close,
  Search,
  ArrowRight,
  Fold,
  Expand,
  Document,
  QuestionFilled,
  Loading
} from '@element-plus/icons-vue'

// Helper to get expand/collapse icon based on section and state
const getExpandIcon = (section: 'stdout' | 'stderr') => {
  return isExpanded[section] ? Fold : Expand
}

// Check if stdout/stderr has actual content (not just whitespace)
const hasContent = (text: string): boolean => {
  return text && text.trim().length > 0
}

// Format output for display
const formatOutput = (text: string): string => {
  return text || '(无输出)'
}

// stderr with tool name prefix for AI analysis context
const stderrWithContext = computed(() => {
  if (!hasContent(props.result.stderr)) return props.result.stderr
  return `[工具: ${props.result.tool_name}]\n${props.result.stderr}`
})

const props = defineProps<{
  result: {
    success: boolean
    tool_name: string
    stdout: string
    stderr: string
    exit_code: number
    file_outputs?: { name: string; filename: string; download_url: string }[]
  }
}>()

const emit = defineEmits<{
  copy: []
  'update:visible': [visible: boolean]
}>()

// Tab state
const activeTab = ref<'summary' | 'files' | 'stdout' | 'stderr'>('summary')
const isMinimized = ref(false)

// Section expansion state
const isExpanded = ref({
  stdout: true,
  stderr: true
})

// Search state
const searchQuery = ref('')
const searchMatches = ref(0)
const searchTotal = ref(0)

// AI analysis state
const aiAnalysis = ref('')
const aiLoading = ref(false)

// Pre references for scrolling
const stdoutPre = ref<HTMLElement | null>(null)
const stderrPre = ref<HTMLElement | null>(null)

// Compute if there's any output
const hasAnyOutput = computed(() => {
  return !!(props.result.stdout || props.result.stderr)
})

// Compute if there are file outputs
const hasFileOutputs = computed(() => {
  return !!(props.result.file_outputs && props.result.file_outputs.length > 0)
})

// Filter STDOUT lines based on search
const filteredStdoutLines = computed(() => {
  if (!searchQuery.value) return props.result.stdout?.split('\n') || []
  const query = searchQuery.value.toLowerCase()
  const lines = props.result.stdout?.split('\n') || []
  searchTotal.value = lines.length

  return lines.filter(line =>
    line.toLowerCase().includes(query)
  )
})

// Filter STDERR lines based on search
const filteredStderrLines = computed(() => {
  if (!searchQuery.value) return props.result.stderr?.split('\n') || []
  const query = searchQuery.value.toLowerCase()
  const lines = props.result.stderr?.split('\n') || []
  searchTotal.value = lines.length

  return lines.filter(line =>
    line.toLowerCase().includes(query)
  )
})

// Switch tab
const switchTab = (tab: 'summary' | 'stdout' | 'stderr') => {
  activeTab.value = tab
  searchQuery.value = ''
  searchMatches.value = 0
}

// Toggle minimize
const toggleMinimize = () => {
  isMinimized.value = !isMinimized.value
  if (!isMinimized.value) {
    nextTick(() => {
      scrollToBottom()
    })
  }
}

// Toggle section expand
const toggleExpand = (section: 'stdout' | 'stderr') => {
  isExpanded.value[section] = !isExpanded.value[section]
}

// Handle search input
const handleSearch = () => {
  searchMatches.value = 0
  if (activeTab.value === 'stdout') {
    searchTotal.value = filteredStdoutLines.value.length
  } else if (activeTab.value === 'stderr') {
    searchTotal.value = filteredStderrLines.value.length
  }
}

// Format output for preview (truncate long lines)
const formatOutputForPreview = (output: string, toolName?: string): string => {
  const prefix = toolName ? `[工具: ${toolName}]\n` : ''
  const lines = (prefix + output).split('\n')
  const maxLines = 10
  const preview = lines.slice(0, maxLines)
  const truncated = lines.length > maxLines
    ? `... (${lines.length - maxLines} more lines)`
    : ''
  return preview.join('\n') + truncated
}

// Copy section
const copySection = async (section: 'stdout' | 'stderr') => {
  let text: string
  if (section === 'stdout') {
    text = props.result.stdout || ''
  } else {
    text = hasContent(props.result.stderr)
      ? `[工具: ${props.result.tool_name}]\n${props.result.stderr}`
      : props.result.stderr || ''
  }
  if (!hasContent(text)) return

  try {
    await navigator.clipboard.writeText(text)
    ElMessage.success(`${section === 'stdout' ? 'STDOUT' : 'STDERR'}已复制到剪贴板`)
  } catch (error) {
    ElMessage.error('复制失败')
  }
}

// Copy all
const copyAll = async () => {
  const text = buildFullResult()
  try {
    await navigator.clipboard.writeText(text)
    ElMessage.success('全部结果已复制到剪贴板')
    emit('copy')
  } catch (error) {
    ElMessage.error('复制失败')
  }
}

// Download section
const downloadSection = (section: 'stdout' | 'stderr', filenameSuffix: string) => {
  const text = section === 'stdout' ? props.result.stdout : props.result.stderr
  if (!hasContent(text)) return

  downloadText(text, `${props.result.tool_name}_${filenameSuffix}.txt`)
}

// Download all
const downloadAll = () => {
  const text = buildFullResult()
  const timestamp = new Date().toISOString().slice(0, 19).replace(/[:T]/g, '-')
  downloadText(text, `${props.result.tool_name}_result_${timestamp}.txt`)
}

// Download file helper
const downloadText = (content: string, filename: string) => {
  const blob = new Blob([content], { type: 'text/plain;charset=utf-8' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
  ElMessage.success(`文件已下载：${filename}`)
}

// Download a single file from server
const downloadFile = async (url: string, filename: string) => {
  try {
    const response = await fetch(url)
    if (!response.ok) {
      ElMessage.error(`下载失败：${filename}`)
      return
    }
    const blob = await response.blob()
    const downloadUrl = URL.createObjectURL(blob)
    const link = document.createElement('a')
    link.href = downloadUrl
    link.download = filename
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
    URL.revokeObjectURL(downloadUrl)
    ElMessage.success(`文件已下载：${filename}`)
  } catch {
    ElMessage.error(`下载失败：${filename}`)
  }
}

// Download all files
const downloadAllFiles = async () => {
  if (!props.result.file_outputs) return
  for (const file of props.result.file_outputs) {
    await downloadFile(file.download_url, file.filename)
  }
}

// Build full result text
const buildFullResult = (): string => {
  const lines = [
    `工具：${props.result.tool_name}`,
    `状态：${props.result.success ? '成功' : '失败'}`,
    `退出码：${props.result.exit_code}`,
    ''
  ]

  if (hasContent(props.result.stdout)) {
    lines.push('--- STDOUT ---')
    lines.push(formatOutput(props.result.stdout))
  }

  if (hasContent(props.result.stderr)) {
    lines.push('')
    lines.push('--- STDERR ---')
    lines.push(formatOutput(props.result.stderr))
  }

  return lines.join('\n')
}

// Scroll to bottom on new result
watch(() => props.result.stdout, (newVal) => {
  if (hasContent(newVal) && !isMinimized.value) {
    nextTick(() => {
      scrollToBottom()
    })
  }
})

watch(() => props.result.stderr, (newVal) => {
  if (hasContent(newVal) && !isMinimized.value) {
    nextTick(() => {
      scrollToBottom()
    })
  }
})

const scrollToBottom = () => {
  // No manual scrolling needed - page uses natural browser scrolling
}

// Ask AI for help analyzing stderr (async polling pattern)
const askAIHelp = async () => {
  if (!props.result.stderr || !hasContent(props.result.stderr)) return

  aiLoading.value = true
  aiAnalysis.value = ''
  try {
    // Step 1: Start the analysis job
    const startRes = await fetch(
      `/api/tools/${props.result.tool_name}/ai-analyze`,
      {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          tool_name: props.result.tool_name,
          stderr: props.result.stderr,
          stdout: props.result.stdout || ''
        })
      }
    )
    const startData = await startRes.json()
    if (!startData.success || !startData.data?.job_id) {
      ElMessage.error(startData.error || '无法启动 AI 分析')
      return
    }

    const jobId = startData.data.job_id

    // Step 2: Poll for the result
    aiLoading.value = false
    aiAnalysis.value = null  // null = showing loading indicator

    const pollInterval = setInterval(async () => {
      const pollRes = await fetch(
        `/api/tools/${props.result.tool_name}/ai-analyze/${jobId}`
      )
      const pollData = await pollRes.json()

      if (!pollData.success) {
        // Job expired or error
        clearInterval(pollInterval)
        aiAnalysis.value = pollData.error || pollData.message || '分析已过期'
        aiAnalysis.value = `[分析失败] ${aiAnalysis.value}`
        ElMessage.warning('AI 分析超时，请稍后重试')
        return
      }

      const status = pollData.data?.status
      if (status === 'running') return  // Still running, keep polling

      // Done or error — stop polling and display result
      clearInterval(pollInterval)

      if (status === 'done') {
        aiAnalysis.value = pollData.data?.result || '(无结果)'
        ElMessage.success('AI 分析完成')
      } else {
        aiAnalysis.value = pollData.data?.result || 'AI 分析失败'
        ElMessage.error('AI 分析失败')
      }
    }, 2000)

    // Step 3: Safety timeout (10 minutes)
    setTimeout(() => {
      clearInterval(pollInterval)
      if (aiAnalysis.value === null) {
        aiAnalysis.value = '[分析超时] AI 分析耗时过长，请稍后查看结果或重试。'
        ElMessage.warning('AI 分析耗时较长，结果可能稍后可用')
      }
    }, 600_000)

  } catch (err) {
    aiAnalysis.value = 'AI 分析请求失败'
    ElMessage.error('AI 分析请求失败')
  }
}

const copyAIAnalysis = async () => {
  if (!aiAnalysis.value) return
  try {
    await navigator.clipboard.writeText(aiAnalysis.value)
    ElMessage.success('已复制到剪贴板')
  } catch {
    ElMessage.error('复制失败')
  }
}

// Render markdown to HTML for AI analysis results
const renderMarkdown = (text: string): string => {
  if (!text) return ''
  // Configure marked for safe rendering
  marked.setOptions({
    breaks: true,
    gfm: true,
  })
  return marked.parse(text) as string
}

// Initialize
nextTick(() => {
  scrollToBottom()
  console.log('=== ResultViewer initialized ===')
  console.log('result.stdout:', props.result.stdout, 'length:', props.result.stdout?.length)
  console.log('result.stderr:', props.result.stderr, 'length:', props.result.stderr?.length)
  console.log('hasContent(stdout):', hasContent(props.result.stdout))
  console.log('hasContent(stderr):', hasContent(props.result.stderr))
})
</script>

<style scoped>
.result-viewer {
  background: #fff;
  padding: 20px;
  border-radius: 8px;
  box-shadow: 0 2px 12px rgba(0, 0, 0, 0.1);
  margin-top: 20px;
  margin-bottom: 20px;
  flex-shrink: 0;
  z-index: 10;
  position: relative;
  scroll-behavior: smooth;
  overflow-y: visible;
  min-height: 100px;
}

.result-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 15px;
  padding-bottom: 12px;
  border-bottom: 2px solid #e4e7ed;
}

.header-left {
  display: flex;
  align-items: center;
  gap: 12px;
}

.result-header h3 {
  margin: 0;
  font-size: 18px;
  color: #303133;
  font-weight: 600;
}

.status-badge {
  font-weight: 500;
}

.exit-code {
  color: #909399;
  font-family: monospace;
}

.header-actions {
  display: flex;
  gap: 8px;
  z-index: 20;
}

.tabs-container {
  margin-bottom: 12px;
}

.tabs-group {
  width: 100%;
}

.tabs-group .el-radio-button__original-radio:checked + .el-radio-button__inner {
  color: #409eff;
  background-color: #ecf5ff;
  border-color: #409eff;
  font-weight: 600;
}

.search-bar {
  margin-bottom: 12px;
  margin-top: 10px;
}

.search-count {
  color: #606266;
  font-size: 12px;
  margin-right: 8px;
  font-weight: 500;
}

.summary-tab {
  padding: 15px 0;
}

.summary-info {
  display: flex;
  flex-direction: column;
  gap: 10px;
  margin-bottom: 20px;
  padding: 15px;
  background: #f5f7fa;
  border-radius: 6px;
}

.info-row {
  display: flex;
  align-items: center;
  gap: 8px;
}

.info-label {
  color: #606266;
  font-weight: 500;
  min-width: 80px;
}

.info-value {
  color: #303133;
  font-weight: 600;
}

.preview-section {
  margin-bottom: 15px;
  padding: 15px;
  border: 1px solid #e4e7ed;
  border-radius: 6px;
  background: #fff;
}

.preview-section.error-section {
  background: #fff;
  border-color: #e4e7ed;
}

.preview-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 10px;
}

.preview-title {
  font-size: 14px;
  font-weight: 600;
  color: #303133;
}

.preview-content {
  font-family: 'Monaco', 'Menlo', monospace;
  font-size: 14px;
  line-height: 1.6;
  color: #303133;
  margin: 0;
  padding: 10px;
  background: #f9f9f9;
  border-radius: 4px;
  white-space: pre;
  word-wrap: normal;
  overflow: auto;
  max-height: 300px;
}

.output-tab {
  display: flex;
  flex-direction: column;
}

.output-controls {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 10px 12px;
}

.control-left, .control-right {
  display: flex;
  gap: 8px;
  align-items: center;
}

.control-left {
  flex: 1;
}

.control-right {
  flex-shrink: 0;
}

.output-title {
  font-size: 14px;
  font-weight: 600;
  color: #303133;
}

.file-list {
  padding: 10px;
}

.file-item {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 8px 12px;
  margin-bottom: 4px;
  background: #f5f7fa;
  border-radius: 4px;
}

.file-name {
  font-family: 'Monaco', 'Menlo', monospace;
  font-size: 13px;
  color: #303133;
}

.output-content {
  padding: 15px;
  background: #fff;
  overflow-x: auto;
  overflow-y: auto;
}

/* Terminal-style styling for code output */
.terminal-output {
  font-family: 'Courier New', 'Courier', monospace;
  background: #1e1e1e;
  color: #d4d4d4;
  padding: 15px;
  border-radius: 4px;
  overflow-x: auto;
  overflow-y: auto;
  white-space: pre;
  word-wrap: normal;
  min-height: 200px;
  display: block;
  line-height: 1.5;
}

.empty-output {
  padding: 40px 20px;
  text-align: center;
}

.preview-close {
  margin-top: 10px;
  text-align: center;
}

.preview-close .el-button {
  color: #909399;
  font-size: 13px;
}

.action-buttons {
  margin-top: 15px;
  padding-top: 12px;
  border-top: 1px solid #e4e7ed;
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
}

.action-buttons .el-button {
  flex: 0 0 auto;
}

/* AI Analysis Section */
.ai-analysis-section {
  margin-top: 15px;
  padding: 15px;
  border: 1px solid #409eff;
  border-radius: 6px;
  background: #ecf5ff;
}

.ai-analysis-section.ai-analysis-loading {
  background: #f0f9ff;
  border-color: #66b1ff;
}

.ai-analysis-status {
  font-size: 12px;
  color: #909399;
  font-style: italic;
}

@keyframes spin {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}

.loading-icon {
  animation: spin 1.5s linear infinite;
  display: inline-block;
}

.ai-analysis-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 10px;
}

.ai-analysis-title {
  font-size: 14px;
  font-weight: 600;
  color: #409eff;
  display: flex;
  align-items: center;
  gap: 4px;
}

.ai-analysis-content {
  background: #fff;
  border-radius: 4px;
  padding: 16px;
  max-height: 500px;
  overflow-y: auto;
  font-size: 14px;
  line-height: 1.8;
  color: #303133;
  word-wrap: break-word;
  white-space: pre-wrap;
}

/* Markdown rendering styles for rendered content */
.ai-analysis-content :deep(h1),
.ai-analysis-content :deep(h2),
.ai-analysis-content :deep(h3),
.ai-analysis-content :deep(h4),
.ai-analysis-content :deep(h5),
.ai-analysis-content :deep(h6) {
  margin: 16px 0 8px 0;
  font-weight: 600;
  color: #1a1a2e;
  line-height: 1.4;
}

.ai-analysis-content :deep(h1) { font-size: 1.4em; }
.ai-analysis-content :deep(h2) { font-size: 1.25em; }
.ai-analysis-content :deep(h3) { font-size: 1.1em; }

.ai-analysis-content :deep(p) {
  margin: 0 0 10px 0;
  line-height: 1.8;
}

.ai-analysis-content :deep(p:last-child) {
  margin-bottom: 0;
}

.ai-analysis-content :deep(ul),
.ai-analysis-content :deep(ol) {
  margin: 8px 0;
  padding-left: 24px;
}

.ai-analysis-content :deep(li) {
  margin-bottom: 4px;
  line-height: 1.7;
}

.ai-analysis-content :deep(li > p) {
  margin: 4px 0;
}

.ai-analysis-content :deep(code) {
  background: #f5f7fa;
  padding: 2px 6px;
  border-radius: 3px;
  font-family: 'Monaco', 'Menlo', monospace;
  font-size: 0.9em;
  color: #e83e8c;
}

.ai-analysis-content :deep(pre) {
  background: #1e1e1e;
  color: #d4d4d4;
  padding: 12px 16px;
  border-radius: 6px;
  overflow-x: auto;
  margin: 10px 0;
  font-family: 'Monaco', 'Menlo', monospace;
  font-size: 13px;
  line-height: 1.6;
}

.ai-analysis-content :deep(pre code) {
  background: transparent;
  padding: 0;
  color: inherit;
}

.ai-analysis-content :deep(blockquote) {
  border-left: 4px solid #409eff;
  margin: 10px 0;
  padding: 8px 16px;
  background: #f0f9ff;
  color: #606266;
}

.ai-analysis-content :deep(blockquote > p) {
  margin: 0;
}

.ai-analysis-content :deep(hr) {
  border: none;
  border-top: 1px solid #e4e7ed;
  margin: 16px 0;
}

.ai-analysis-content :deep(strong) {
  font-weight: 600;
  color: #1a1a2e;
}

.ai-analysis-content :deep(a) {
  color: #409eff;
  text-decoration: none;
}

.ai-analysis-content :deep(a:hover) {
  text-decoration: underline;
}

.ai-analysis-content :deep(table) {
  border-collapse: collapse;
  width: 100%;
  margin: 10px 0;
  font-size: 13px;
}

.ai-analysis-content :deep(th),
.ai-analysis-content :deep(td) {
  border: 1px solid #e4e7ed;
  padding: 8px 12px;
  text-align: left;
}

.ai-analysis-content :deep(th) {
  background: #f5f7fa;
  font-weight: 600;
}

.preview-header-actions {
  display: flex;
  align-items: center;
  gap: 8px;
}
</style>
