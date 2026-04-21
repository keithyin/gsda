<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>Low Q Analysis - 低质量产率问题分析</h3>
      <p class="description">分析低 Q20/Q30 产率问题，综合执行多项分析：subreads-SMC 一致性检查、正反链不一致检测、homo+STR 区域覆盖率和Macebell reads比例</p>
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
          placeholder="格式：user@host"
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
    </div>

    <!-- SBR BAM File -->
    <div class="form-group">
      <label>SBR BAM 文件 <span class="required">*</span></label>
      <div v-if="fileSource === 'local'" class="local-file-section">
        <!-- Local File Paths -->
        <div class="input-row">
          <el-input
            v-model="formData.sbr_bam"
            placeholder="/path/to/sbr.bam"
            clearable
          ></el-input>
        </div>

        <!-- File Upload -->
        <div
          class="file-upload-zone"
          @dragover="handleDragOver"
          @dragleave="handleDragLeave"
          @drop="handleDrop('sbr_bam')"
          @click="triggerUpload('sbr')"
        >
          <div v-if="uploadedFiles.sbr.length === 0">
            <el-icon :size="48" color="#909399"><Upload /></el-icon>
            <p style="margin-top: 10px; color: #909399;">拖拽 BAM 文件到此处，或点击上传</p>
          </div>
          <div v-else>
            <el-tag
              v-for="(file, index) in uploadedFiles.sbr"
              :key="'sbr-' + index"
              closable
              @close="removeUploadedFile('sbr', index)"
              type="success"
              style="margin-right: 8px; margin-bottom: 8px;"
            >
              {{ file.name }}
            </el-tag>
            <el-button size="small" type="primary" @click="triggerUpload('sbr')">
              <el-icon><Plus /></el-icon> 添加文件
            </el-button>
          </div>
          <input
            ref="sbrUploadInput"
            type="file"
            accept=".bam"
            style="display: none"
            @change="handleUpload('sbr')"
          />
        </div>
      </div>

      <!-- Remote File Input -->
      <div v-else class="remote-file-section">
        <el-input
          v-model="formData.sbr_bam"
          placeholder="user@host:/path/to/sbr.bam"
          clearable
        ></el-input>
      </div>
    </div>

    <!-- SMC BAM File -->
    <div class="form-group">
      <label>SMC BAM 文件 <span class="required">*</span></label>
      <div v-if="fileSource === 'local'" class="local-file-section">
        <!-- Local File Paths -->
        <div class="input-row">
          <el-input
            v-model="formData.smc_bam"
            placeholder="/path/to/smc.bam"
            clearable
          ></el-input>
        </div>

        <!-- File Upload -->
        <div
          class="file-upload-zone"
          @dragover="handleDragOver"
          @dragleave="handleDragLeave"
          @drop="handleDrop('smc_bam')"
          @click="triggerUpload('smc')"
        >
          <div v-if="uploadedFiles.smc.length === 0">
            <el-icon :size="48" color="#909399"><Upload /></el-icon>
            <p style="margin-top: 10px; color: #909399;">拖拽 BAM 文件到此处，或点击上传</p>
          </div>
          <div v-else>
            <el-tag
              v-for="(file, index) in uploadedFiles.smc"
              :key="'smc-' + index"
              closable
              @close="removeUploadedFile('smc', index)"
              type="success"
              style="margin-right: 8px; margin-bottom: 8px;"
            >
              {{ file.name }}
            </el-tag>
            <el-button size="small" type="primary" @click="triggerUpload('smc')">
              <el-icon><Plus /></el-icon> 添加文件
            </el-button>
          </div>
          <input
            ref="smcUploadInput"
            type="file"
            accept=".bam"
            style="display: none"
            @change="handleUpload('smc')"
          />
        </div>
      </div>

      <!-- Remote File Input -->
      <div v-else class="remote-file-section">
        <el-input
          v-model="formData.smc_bam"
          placeholder="user@host:/path/to/smc.bam"
          clearable
        ></el-input>
      </div>
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

const fileSource = ref<'local' | 'remote'>('local')
const sshConfig = ref({
  server: '',
  password: ''
})
const uploadedFiles = reactive({
  sbr: [] as {name: string, localPath: string}[],
  smc: [] as {name: string, localPath: string}[]
})
const sbrUploadInput = ref<HTMLInputElement | null>(null)
const smcUploadInput = ref<HTMLInputElement | null>(null)

const formData = reactive({
  sbr_bam: '',
  smc_bam: ''
})

const loading = ref(false)

// Sync loading state with parent's isExecuting prop
watch(() => props.isExecuting, (newVal) => {
  loading.value = newVal
})

// File upload handlers for SBR
const triggerUpload = (type: 'sbr' | 'smc') => {
  if (type === 'sbr') {
    sbrUploadInput.value?.click()
  } else {
    smcUploadInput.value?.click()
  }
}

const handleUpload = async (type: 'sbr' | 'smc') => {
  const inputRef = type === 'sbr' ? sbrUploadInput : smcUploadInput
  const target = inputRef.value as HTMLInputElement
  const files = target?.files
  if (!files) return

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    if (!file.name.endsWith('.bam')) {
      ElMessage.warning(`只支持 .bam 文件：${file.name}`)
      continue
    }

    try {
      const formDataObj = new FormData()
      formDataObj.append('file', file)

      const response = await axios.post('/api/files/upload', formDataObj)
      if (response.data.success) {
        if (type === 'sbr') {
          uploadedFiles.sbr.push({
            name: file.name,
            localPath: response.data.data.local_path
          })
          formData.sbr_bam = response.data.data.local_path
        } else {
          uploadedFiles.smc.push({
            name: file.name,
            localPath: response.data.data.local_path
          })
          formData.smc_bam = response.data.data.local_path
        }
      }
    } catch (error: any) {
      ElMessage.error(`上传失败：${error.message || '未知错误'}`)
    }
  }

  // Reset input
  if (inputRef.value) {
    inputRef.value.value = ''
  }
}

const removeUploadedFile = (type: 'sbr' | 'smc', index: number) => {
  const file = (uploadedFiles as any)[type][index]
  (uploadedFiles as any)[type].splice(index, 1)
  // Clear the field if we removed the only file
  if ((uploadedFiles as any)[type].length === 0) {
    if (type === 'sbr') {
      formData.sbr_bam = ''
    } else {
      formData.smc_bam = ''
    }
  }
}

// File upload handlers for SMC

// Drag and drop handlers
const handleDragOver = (e: DragEvent) => {
  e.preventDefault()
}

const handleDragLeave = (e: DragEvent) => {
  e.preventDefault()
}

const handleDrop = (type: 'sbr_bam' | 'smc_bam') => async (e: DragEvent) => {
  e.preventDefault()
  const files = e.dataTransfer?.files
  if (!files) return

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    if (!file.name.endsWith('.bam')) {
      ElMessage.warning(`只支持 .bam 文件：${file.name}`)
      continue
    }

    try {
      const formDataObj = new FormData()
      formDataObj.append('file', file)

      const response = await axios.post('/api/files/upload', formDataObj)
      if (response.data.success) {
        if (type === 'sbr_bam') {
          uploadedFiles.sbr.push({
            name: file.name,
            localPath: response.data.data.local_path
          })
          formData.sbr_bam = response.data.data.local_path
        } else {
          uploadedFiles.smc.push({
            name: file.name,
            localPath: response.data.data.local_path
          })
          formData.smc_bam = response.data.data.local_path
        }
      }
    } catch (error: any) {
      ElMessage.error(`上传失败：${error.message || '未知错误'}`)
    }
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
  console.log('=== LowQAnalysisForm - execute() START ===')
  console.log('formData:', formData)
  console.log('fileSource:', fileSource.value)
  console.log('sshConfig:', sshConfig.value)

  // Validate required fields
  if (!formData.sbr_bam) {
    ElMessage.error('请选择或输入 SBR BAM 文件')
    return
  }
  if (!formData.smc_bam) {
    ElMessage.error('请选择或输入 SMC BAM 文件')
    return
  }

  // Scroll to ensure button is visible
  await nextTick()
  scrollToActions()

  loading.value = true
  try {
    // Build request body
    const request: any = {
      tool_name: 'low-q-analysis',
      sbr_bam: formData.sbr_bam,
      smc_bam: formData.smc_bam
    }

    // Add SSH config for remote files
    if (fileSource.value === 'remote') {
      request.ssh_server = sshConfig.value.server
      request.ssh_password = sshConfig.value.password
    }

    console.log('=== Emitting execute event ===')
    console.log('Full request object:', JSON.stringify(request, null, 2))

    emit('execute', request)
  } catch (error) {
    console.error('=== execute() ERROR ===', error)
    throw error
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

.required {
  color: #f56c6c;
  margin-left: 4px;
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

.local-file-section {
  margin-bottom: 15px;
}

.remote-file-section {
  margin-bottom: 15px;
}

.form-actions {
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e4e7ed;
  flex-shrink: 0;
}
</style>
