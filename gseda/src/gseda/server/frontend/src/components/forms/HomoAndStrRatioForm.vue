<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>HomoAndStr Ratio - 同源与串联重复比率分析</h3>
      <p class="description">统计 BAM/FASTQ/FASTA 文件中同源重复（homo）和串联重复（str）区域的比率。</p>
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
        <label>文件路径 (远程)</label>
        <div v-for="(path, index) in formData.files" :key="'remote-' + index" class="input-row">
          <el-input
            v-model="formData.files[index]"
            placeholder="user@host:/path/to/file.bam 或 /path/to/*.fq"
            clearable
          >
            <template #append>
              <el-button @click="removeFile(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addFile">
          <el-icon><Plus /></el-icon> 添加文件
        </el-button>
      </div>
    </div>

    <!-- Local File Inputs -->
    <div v-else>
      <!-- File Paths -->
      <div class="form-group">
        <label>文件路径 (本地)</label>
        <div v-for="(path, index) in formData.files" :key="'local-' + index" class="input-row">
          <el-input
            v-model="formData.files[index]"
            placeholder="/path/to/file.bam 或 /path/to/*.fq"
            clearable
          >
            <template #append>
              <el-button @click="removeFile(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addFile">
          <el-icon><Plus /></el-icon> 添加文件
        </el-button>
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
            <p style="margin-top: 10px; color: #909399;">拖拽 .bam/.fq/.fastq/.fa/.fasta 文件到此处，或点击上传</p>
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
            accept=".bam,.fq,.fastq,.fa,.fasta"
            style="display: none"
            @change="handleUpload"
            multiple
          />
        </div>
      </div>
    </div>

    <!-- RQ Threshold -->
    <div class="form-group">
      <label>Minimum RQ (可选)</label>
      <el-input-number
        v-model="formData.rq_thr"
        :min="0"
        :max="1"
        :step="0.05"
        :precision="2"
        :placeholder="0.95"
      >
        <template #append>
          <span>RQ</span>
        </template>
      </el-input-number>
      <p class="hint">过滤 reads 的最小质量值，默认 0.95</p>
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
const uploadedFiles = ref<{name: string, localPath: string}[]>([])
const uploadInput = ref<HTMLInputElement | null>(null)

const formData = reactive({
  files: ['/path/to/file1.bam'],
  rq_thr: 0.95
})

const loading = ref(false)

// Sync loading state with parent's isExecuting prop
watch(() => props.isExecuting, (newVal) => {
  loading.value = newVal
})

const addFile = () => {
  formData.files.push('/path/to/file.bam')
}

const removeFile = (index: number) => {
  if (formData.files.length > 1) {
    formData.files.splice(index, 1)
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

  const validExtensions = ['.bam', '.fq', '.fastq', '.fa', '.fasta']

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    const hasValidExt = validExtensions.some(ext => file.name.toLowerCase().endsWith(ext))

    if (!hasValidExt) {
      ElMessage.warning(`只支持 .bam, .fq, .fastq, .fa, .fasta 文件: ${file.name}`)
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
        formData.files.push(response.data.data.local_path)
      }
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
    }
  }

  // Reset input
  if (uploadInput.value) {
    uploadInput.value.value = ''
  }
}

const removeUploadedFile = (index: number) => {
  const file = uploadedFiles.value[index]
  uploadedFiles.value.splice(index, 1)
  // Remove from files array
  formData.files = formData.files.filter(p => p !== file.localPath)
}

// Drag and drop handlers
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

  const validExtensions = ['.bam', '.fq', '.fastq', '.fa', '.fasta']

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    const hasValidExt = validExtensions.some(ext => file.name.toLowerCase().endsWith(ext))

    if (!hasValidExt) {
      ElMessage.warning(`只支持 .bam, .fq, .fastq, .fa, .fasta 文件: ${file.name}`)
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
        formData.files.push(response.data.data.local_path)
      }
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
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
  console.log('=== HomoAndStrRatioForm - execute() START ===')
  console.log('formData:', formData)
  console.log('fileSource:', fileSource.value)
  console.log('sshConfig:', sshConfig.value)

  // Scroll to ensure button is visible
  await nextTick()
  scrollToActions()

  loading.value = true
  try {
    // Build request body
    const request: any = {
      tool_name: 'homo-and-str-ratio',
      files: formData.files,
      rq_thr: formData.rq_thr
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
