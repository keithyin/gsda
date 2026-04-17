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
            <el-button link @click="switchTab('stderr')">
              查看完整 <el-icon><ArrowRight /></el-icon>
            </el-button>
          </div>
          <pre class="preview-content">{{ formatOutputForPreview(result.stderr) }}</pre>
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
        <div class="output-content error-bg terminal-output">
          <pre ref="stderrPre" class="output-text">
            {{ result.stderr }}
          </pre>
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
import {
  DocumentCopy,
  Download,
  Close,
  Search,
  ArrowRight,
  Fold,
  Expand,
  Document
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

const props = defineProps<{
  result: {
    success: boolean
    tool_name: string
    stdout: string
    stderr: string
    exit_code: number
  }
}>()

const emit = defineEmits<{
  copy: []
  'update:visible': [visible: boolean]
}>()

// Tab state
const activeTab = ref<'summary' | 'stdout' | 'stderr'>('summary')
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

// Pre references for scrolling
const stdoutPre = ref<HTMLElement | null>(null)
const stderrPre = ref<HTMLElement | null>(null)

// Compute if there's any output
const hasAnyOutput = computed(() => {
  return !!(props.result.stdout || props.result.stderr)
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
const formatOutputForPreview = (output: string): string => {
  const lines = output.split('\n')
  const maxLines = 10
  const preview = lines.slice(0, maxLines)
  const truncated = lines.length > maxLines
    ? `... (${lines.length - maxLines} more lines)`
    : ''
  return preview.join('\n') + truncated
}

// Copy section
const copySection = async (section: 'stdout' | 'stderr') => {
  const text = section === 'stdout' ? props.result.stdout : props.result.stderr
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
  background: #fef0f0;
  border-color: #fde2e2;
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

.output-content {
  padding: 15px;
  background: #fff;
  overflow-x: auto;
  overflow-y: auto;
}

.output-content.error-bg {
  background: #fef0f0;
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
</style>
