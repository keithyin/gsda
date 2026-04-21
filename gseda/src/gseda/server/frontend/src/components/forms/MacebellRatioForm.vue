<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>Macebell Ratio - Macebell 比例分析</h3>
      <p class="description">检测 Macebell reads 比例</p>
    </div>

    <!-- File Source Toggle -->
    <div class="form-group">
      <el-radio-group v-model="fileSource" size="medium">
        <el-radio-button label="local">本地文件</el-radio-button>
        <el-radio-button label="remote">远程文件 (SCP)</el-radio-button>
      </el-radio-group>
    </div>

    <!-- Remote File Inputs -->
    <div v-if="fileSource === 'remote'" class="remote-config">
      <div class="form-group">
        <label>SSH 服务器地址</label>
        <el-input
          v-model="sshConfig.server"
          placeholder="格式: user@host"
          clearable
        >
          <template #prepend>SSH</template>
        </el-input>
      </div>

      <div class="form-group">
        <label>SSH 密码</label>
        <el-input
          v-model="sshConfig.password"
          type="password"
          placeholder="SSH 密码"
          show-password
        ></el-input>
      </div>

      <div class="form-group">
        <label>BAM 文件路径 (远程)</label>
        <el-input
          v-model="formData.input_bam"
          placeholder="user@host:/path/to/file.bam"
          clearable
        >
          <template #prepend>文件</template>
        </el-input>
      </div>
    </div>

    <!-- Local File Inputs -->
    <div v-else>
      <!-- File Path -->
      <div class="form-group">
        <label>BAM 文件路径 (本地)</label>
        <el-input
          v-model="formData.input_bam"
          placeholder="/path/to/file.bam"
          clearable
        >
          <template #prepend>文件</template>
        </el-input>
      </div>

      <!-- File Upload -->
      <div class="form-group">
        <label>或上传文件</label>
        <div
          class="file-upload-zone"
          @dragover="handleDragOver"
          @dragleave="handleDragLeave"
          @drop="handleDrop"
          @click="triggerUpload"
        >
          <div v-if="uploadedFiles.length === 0">
            <el-icon :size="48" color="#909399"><Upload /></el-icon>
            <p style="margin-top: 10px; color: #909399;">拖拽 .bam 文件到此处，或点击上传</p>
          </div>
          <div v-else>
            <el-tag
              v-for="(file, index) in uploadedFiles"
              :key="index"
              closable
              @close="removeUploadedFile(index)"
              type="success"
              style="margin-right: 8px; margin-bottom: 8px;"
            >
              {{ file.name }}
            </el-tag>
            <el-button size="small" type="primary" @click="triggerUpload">
              <el-icon><Plus /></el-icon> 添加文件
            </el-button>
          </div>
          <input
            ref="uploadInput"
            type="file"
            accept=".bam"
            style="display: none"
            @change="handleUpload"
          />
        </div>
      </div>
    </div>

    <!-- Threads -->
    <div class="form-group">
      <label>线程数 (可选)</label>
      <el-input-number
        v-model="formData.threads"
        :min="1"
        :step="1"
        :placeholder="16"
      >
        <template #append>
          <span>threads</span>
        </template>
      </el-input-number>
      <p class="hint">消费者进程数，默认使用所有 CPU 核心</p>
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
import { Plus, Upload } from '@element-plus/icons-vue'
import axios from 'axios'
import { ElMessage } from 'element-plus'

const emit = defineEmits<{
  execute: [data: any]
}>()

const props = defineProps<{
  isExecuting?: boolean
}>()

const fileSource = ref<'local' | 'remote'>('local')
const sshConfig = ref({
  server: '',
  password: ''
})
const uploadedFiles = ref<{name: string, localPath: string}[]>([])
const uploadInput = ref<HTMLInputElement | null>(null)

const formData = reactive({
  input_bam: '/path/to/file.bam',
  threads: 16
})

const loading = ref(false)

// Sync loading state with parent's isExecuting prop
watch(() => props.isExecuting, (newVal) => {
  loading.value = newVal
})

const triggerUpload = () => {
  uploadInput.value?.click()
}

const handleUpload = async (event: Event) => {
  const target = event.target as HTMLInputElement
  const files = target.files
  if (!files) return

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    if (!file.name.toLowerCase().endsWith('.bam')) {
      ElMessage.warning('只支持 .bam 文件: ' + file.name)
      continue
    }

    try {
      const formDataObj = new FormData()
      formDataObj.append('file', file)

      const response = await axios.post('/api/files/upload', formDataObj)
      if (response.data.success) {
        uploadedFiles.value.push({
          name: file.name,
          localPath: response.data.data.local_path
        })
        formData.input_bam = response.data.data.local_path
      }
    } catch (error: any) {
      ElMessage.error('上传失败: ' + (error.message || '未知错误'))
    }
  }

  if (uploadInput.value) {
    uploadInput.value.value = ''
  }
}

const removeUploadedFile = (index: number) => {
  const file = uploadedFiles.value[index]
  uploadedFiles.value.splice(index, 1)
  formData.input_bam = '/path/to/file.bam'
}

const handleDragOver = (e: DragEvent) => {
  e.preventDefault()
}

const handleDragLeave = (e: DragEvent) => {
  e.preventDefault()
}

const handleDrop = async (e: DragEvent) => {
  e.preventDefault()
  const files = e.dataTransfer?.files
  if (!files) return

  const file = files[0]
  if (!file.name.toLowerCase().endsWith('.bam')) {
    ElMessage.warning('只支持 .bam 文件')
    return
  }

  try {
    const formDataObj = new FormData()
    formDataObj.append('file', file)

    const response = await axios.post('/api/files/upload', formDataObj)
    if (response.data.success) {
      uploadedFiles.value.push({
        name: file.name,
        localPath: response.data.data.local_path
      })
      formData.input_bam = response.data.data.local_path
    }
  } catch (error: any) {
    ElMessage.error('上传失败: ' + (error.message || '未知错误'))
  }
}

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
      tool_name: 'macebell-ratio',
      input_bam: formData.input_bam,
      threads: formData.threads
    }

    if (fileSource.value === 'remote') {
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

.file-upload-zone.drag-over {
  border-color: #409eff;
  background: #ecf5ff;
}

.remote-config {
  background: #f0f9ff;
  padding: 20px;
  border-radius: 6px;
  border: 1px solid #bae6fd;
  margin-bottom: 20px;
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
</style>
