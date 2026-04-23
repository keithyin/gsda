<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>BAM Basic Stat - BAM 基本质量统计</h3>
      <p class="description">统计 BAM 文件中 reads 的基本质量指标，包括 channel 分布、读长统计、质量分数分布（Q8~Q30）等</p>
    </div>

    <!-- BAM File - Unified Section -->
    <div class="form-group bam-unity-section">
      <label>BAM 文件来源</label>
      <el-radio-group v-model="fileSource" size="medium">
        <el-radio-button value="local">服务器本地</el-radio-button>
        <el-radio-button value="upload">客户端上传</el-radio-button>
        <el-radio-button value="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>

      <!-- Server Local -->
      <div v-if="fileSource === 'local'" class="bam-input-area">
        <label class="sub-label">BAM 文件路径 (支持通配符*)</label>
        <div v-for="(path, index) in formData.bams" :key="'local-' + index" class="input-row">
          <el-input v-model="formData.bams[index]"
            placeholder="/path/to/*.bam 或 /path/to/file.bam" clearable>
            <template #append>
              <el-button @click="removeBam(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addBam">
          <el-icon><Plus /></el-icon> 添加 BAM 文件
        </el-button>
      </div>

      <!-- Client Upload -->
      <div v-else-if="fileSource === 'upload'" class="bam-input-area">
        <div class="file-upload-zone" @click="triggerUpload">
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 BAM 文件</p>
        </div>
        <input ref="uploadInput" type="file" accept=".bam"
          style="display: none" @change="handleUpload" />
      </div>

      <!-- SCP Remote -->
      <div v-else-if="fileSource === 'scp'" class="bam-input-area bam-scp-area">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input v-model="sshConfig.server" placeholder="格式: user@host" clearable
            autocomplete="url" />
        </div>
        <div class="form-group">
          <label>SSH 密码</label>
          <el-input v-model="sshConfig.password" type="password" placeholder="SSH 密码" show-password
            autocomplete="current-password" />
        </div>
        <div class="form-group">
          <label>BAM 文件路径 (支持通配符*)</label>
          <div v-for="(path, index) in formData.bams" :key="'remote-' + index" class="input-row">
            <el-input v-model="formData.bams[index]"
              placeholder="/path/to/*.bam 或 /path/to/file.bam" clearable>
              <template #append>
                <el-button @click="removeBam(index)" :disabled="index === 0">
                  <el-icon><Minus /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>
          <el-button type="primary" plain @click="addBam">
            <el-icon><Plus /></el-icon> 添加 BAM 文件
          </el-button>
        </div>
      </div>
    </div>

    <!-- Channel Tag -->
    <div class="form-group">
      <label>Channel Tag</label>
      <el-select v-model="formData.channel_tag" placeholder="选择 channel tag">
        <el-option label="ch (标准 channel)" value="ch" />
        <el-option label="zm (Zymo 标准)" value="zm" />
      </el-select>
    </div>

    <!-- Min RQ -->
    <div class="form-group">
      <label>Minimum RQ (可选)</label>
      <el-input-number
        v-model="formData.min_rq"
        :min="0"
        :max="1"
        :step="0.1"
        placeholder="留空则不设置最低质量标准"
      >
        <template #append>
          <span>RQ</span>
        </template>
      </el-input-number>
    </div>

    <!-- Execute Button -->
    <div class="form-actions">
      <el-button
        type="primary"
        size="large"
        :loading="loading"
        @click="execute"
      >
        {{ loading ? '分析中...' : '开始分析' }}
      </el-button>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, reactive, watch, nextTick } from 'vue'
import { Plus, Minus, Upload } from '@element-plus/icons-vue'
import axios from 'axios'
import { ElMessage } from 'element-plus'

const emit = defineEmits<{
  execute: [data: any]
}>()

const props = defineProps<{
  isExecuting?: boolean
}>()

const fileSource = ref<'local' | 'upload' | 'scp'>('local')
const sshConfig = ref({
  server: '',
  password: ''
})
const uploadInput = ref<HTMLInputElement | null>(null)

const formData = reactive({
  bams: ['/path/to/file1.bam'],
  channel_tag: 'ch',
  min_rq: null as number | null
})

const loading = ref(false)

// Sync loading state with parent's isExecuting prop
watch(() => props.isExecuting, (newVal) => {
  loading.value = newVal
})

const addBam = () => {
  formData.bams.push('/path/to/file.bam')
}

const removeBam = (index: number) => {
  if (formData.bams.length > 1) {
    formData.bams.splice(index, 1)
  }
}

// File upload handlers
const triggerUpload = () => {
  uploadInput.value?.click()
}

const handleUpload = async (event: Event) => {
  const target = event.target as HTMLInputElement
  const files = target.files
  if (!files) return

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    if (!file.name.endsWith('.bam')) {
      ElMessage.warning(`只支持 .bam 文件: ${file.name}`)
      continue
    }

    try {
      const formDataObj = new FormData()
      formDataObj.append('file', file)

      const response = await axios.post('/api/files/upload', formDataObj)
      if (response.data.success) {
        formData.bams.push(response.data.data.local_path)
      }
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
    }
  }

  if (uploadInput.value) {
    uploadInput.value.value = ''
  }
}

// Helper to scroll to form actions
const scrollToActions = () => {
  const formContainer = document.querySelector('.form-container') as HTMLElement
  if (formContainer) {
    const actions = formContainer.querySelector('.form-actions') as HTMLElement
    if (actions) {
      actions.scrollIntoView({ behavior: 'smooth', block: 'center' })
    }
  }
}

const execute = async () => {
  await nextTick()
  scrollToActions()

  loading.value = true
  try {
    const request: any = {
      tool_name: 'bam-basic-stat',
      bams: formData.bams,
      channel_tag: formData.channel_tag,
      min_rq: formData.min_rq
    }

    if (fileSource.value === 'scp') {
      request.ssh_server = sshConfig.value.server
      request.ssh_password = sshConfig.value.password
    }

    emit('execute', request)
  } finally {
    loading.value = false
  }
}
</script>

<style scoped>
.form-container {
  background: #fff;
  padding: 25px;
  border-radius: 8px;
  box-shadow: 0 2px 12px rgba(0, 0, 0, 0.1);
  overflow-y: auto;
  display: flex;
  flex-direction: column;
  transition: all 0.3s;
  min-height: 0;
  flex: 1;
  max-height: calc(100vh - 400px);
  z-index: 1;
}

.form-container.is-loading {
  opacity: 0.6;
  pointer-events: none;
  position: relative;
}

.form-container.is-loading::after {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(255, 255, 255, 0.5);
  display: flex;
  justify-content: center;
  align-items: center;
  font-size: 18px;
  color: #409eff;
  font-weight: 500;
}

h3 {
  margin: 0 0 10px 0;
  color: #303133;
  font-size: 18px;
}

.description {
  color: #909399;
  font-size: 14px;
  margin-bottom: 25px;
}

.form-section {
  border-bottom: 1px solid #e4e7ed;
  padding-bottom: 20px;
  margin-bottom: 25px;
}

.form-group {
  margin-bottom: 25px;
}

.form-group label {
  display: block;
  margin-bottom: 10px;
  color: #606266;
  font-weight: 500;
}

.input-row {
  margin-bottom: 10px;
}

.file-upload-zone {
  border: 2px dashed #d9d9d9;
  border-radius: 6px;
  padding: 40px 20px;
  text-align: center;
  cursor: pointer;
  transition: all 0.3s;
  position: relative;
}

.file-upload-zone:hover {
  border-color: #409eff;
  background: #f5f7fa;
}

.form-actions {
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e4e7ed;
  flex-shrink: 0;
}

.hint {
  font-size: 12px;
  color: #909399;
  margin-top: 6px;
  margin-left: 1px;
}

.bam-unity-section {
  border: 1px solid #d9d9d9;
  border-radius: 8px;
  padding: 20px;
  background: #fafbfc;
}

.bam-unity-section .el-radio-group {
  margin-bottom: 16px;
}

.bam-input-area {
  margin-top: 16px;
  padding-top: 16px;
  border-top: 1px dashed #e4e7ed;
}

.bam-input-area .sub-label {
  display: block;
  margin-bottom: 10px;
  color: #606266;
  font-weight: 500;
}

.bam-scp-area {
  background: #f0f9ff;
  padding: 16px 20px;
  border-radius: 6px;
  border: 1px solid #bae6fd;
}
</style>
