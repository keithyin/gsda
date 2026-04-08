<template>
  <div class="result-viewer">
    <!-- Result Header -->
    <div class="result-header">
      <h3>执行结果</h3>
      <div class="header-actions">
        <el-tag :type="result.success ? 'success' : 'danger'" size="small">
          {{ result.success ? '成功' : '失败' }}
        </el-tag>
        <span class="exit-code" style="margin-left: 10px; font-size: 13px;">
          退出码：{{ result.exit_code }}
        </span>
      </div>
    </div>

    <!-- Result Content -->
    <div class="result-output">
      <div v-if="result.stdout" class="output-section">
        <div class="section-title">STDOUT</div>
        <pre v-html="result.stdout"></pre>
      </div>

      <div v-if="result.stderr" class="output-section error-section">
        <div class="section-title">STDERR (错误信息)</div>
        <pre v-html="result.stderr"></pre>
      </div>

      <div v-if="!result.stdout && !result.stderr" class="empty-output">
        无输出内容
      </div>
    </div>

    <!-- Action Buttons -->
    <div class="action-buttons">
      <el-button @click="copyResult">
        <el-icon><DocumentCopy /></el-icon>
        复制结果
      </el-button>
      <el-button @click="$emit('update:visible', false)">关闭</el-button>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue'
import { ElMessage } from 'element-plus'
import { DocumentCopy } from '@element-plus/icons-vue'

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

// Build full result text for copying
const fullResult = computed(() => {
  const lines = [
    `工具：${props.result.tool_name}`,
    `状态：${props.result.success ? '成功' : '失败'}`,
    `退出码：${props.result.exit_code}`,
    ''
  ]

  if (props.result.stdout) {
    lines.push('--- STDOUT ---')
    lines.push(props.result.stdout)
  }

  if (props.result.stderr) {
    lines.push('')
    lines.push('--- STDERR ---')
    lines.push(props.result.stderr)
  }

  return lines.join('\n')
})

// Copy result to clipboard
const copyResult = () => {
  navigator.clipboard.writeText(fullResult.value).then(() => {
    ElMessage.success('结果已复制到剪贴板')
    emit('copy')
  })
}
</script>

<style scoped>
.result-viewer {
  background: #fff;
  padding: 25px;
  border-radius: 8px;
  box-shadow: 0 2px 12px rgba(0, 0, 0, 0.1);
  margin-top: 30px;
}

.result-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
  padding-bottom: 15px;
  border-bottom: 1px solid #e4e7ed;
}

.result-header h3 {
  margin: 0;
  font-size: 18px;
  color: #303133;
}

.header-actions {
  display: flex;
  align-items: center;
}

.exit-code {
  color: #909399;
}

.result-output {
  background: #282c34;
  border-radius: 6px;
  max-height: 500px;
  overflow-y: auto;
}

.output-section {
  margin-bottom: 20px;
  padding: 15px;
}

.output-section:last-child {
  margin-bottom: 0;
}

.section-title {
  color: #61affe;
  font-size: 13px;
  font-weight: 600;
  margin-bottom: 10px;
  text-transform: uppercase;
}

.output-section pre {
  color: #abb2bf;
  font-family: 'Monaco', 'Menlo', monospace;
  font-size: 13px;
  line-height: 1.5;
  margin: 0;
  white-space: pre-wrap;
  word-wrap: break-word;
}

.error-section {
  background: rgba(245, 108, 108, 0.1);
}

.error-section .section-title {
  color: #f56c6c;
}

.empty-output {
  padding: 30px;
  text-align: center;
  color: #909399;
}

.action-buttons {
  margin-top: 20px;
  padding-top: 15px;
  border-top: 1px solid #e4e7ed;
  display: flex;
  justify-content: flex-end;
  gap: 10px;
}
</style>
