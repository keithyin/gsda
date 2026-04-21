<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>Reads Quality HP </h3>
      <p class="description">使用 gsmm2-metric (hp-v2) 分析 BAM 中 poly-N 区域的质量</p>
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
        <label>BAM 文件路径 (支持通配符*)</label>
        <div v-for="(path, index) in formData.bams" :key="'remote-' + index" class="input-row">
          <el-input
            v-model="formData.bams[index]"
            placeholder="/path/to/*.bam 或 /path/to/file.bam"
            clearable
          >
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

    <!-- Local File Inputs -->
    <div v-else>
      <!-- BAM File Paths -->
      <div class="form-group">
        <label>BAM 文件路径 (支持通配符*)</label>
        <div v-for="(path, index) in formData.bams" :key="'local-' + index" class="input-row">
          <el-input
            v-model="formData.bams[index]"
            placeholder="/path/to/*.bam 或 /path/to/file.bam"
            clearable
          >
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

      <!-- File Upload -->
      <div class="form-group">
        <label>或上传单个 BAM 文件</label>
        <div
          class="file-upload-zone"
          @click="triggerUpload"
        >
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 BAM 文件</p>
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

    <!-- Reference FASTA Files -->
    <div class="form-group">
      <label>Reference FASTA 文件 (必填)</label>
      <div v-for="(path, index) in formData.refs" :key="'ref-' + index" class="input-row">
        <el-input
          v-model="formData.refs[index]"
          placeholder="/path/to/reference.fa"
          clearable
        >
          <template #append>
            <el-button @click="removeRef(index)" :disabled="index === 0">
              <el-icon><Minus /></el-icon>
            </el-button>
          </template>
        </el-input>
      </div>
      <el-button type="primary" plain @click="addRef">
        <el-icon><Plus /></el-icon> 添加 Reference
      </el-button>
      <p class="hint">如果不提供，则无法生成对齐相关的指标。单个 Reference 可自动对应多个 BAM 文件</p>
    </div>

    <!-- Reference FASTA Source Toggle -->
    <div class="form-group">
      <label>Reference FASTA 文件来源</label>
      <el-radio-group v-model="refSource" size="medium">
        <el-radio-button label="local">服务器本地</el-radio-button>
        <el-radio-button label="upload">客户端上传</el-radio-button>
        <el-radio-button label="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>
    </div>

    <!-- Remote Reference Files (SCP) -->
    <div v-if="refSource === 'scp'" class="remote-config">
      <div class="form-group">
        <label>SSH 服务器地址</label>
        <el-input
          v-model="scpRefConfig.server"
          placeholder="格式：user@host"
          clearable
        >
          <template #prepend>SSH</template>
        </el-input>
      </div>

      <div class="form-group">
        <label>SSH 密码</label>
        <el-input
          v-model="scpRefConfig.password"
          type="password"
          placeholder="SSH 密码"
          show-password
        ></el-input>
      </div>

      <div class="form-group">
        <label>Reference FASTA 文件路径 (支持通配符*)</label>
        <div v-for="(path, index) in formData.scp_refs" :key="'scp-ref-' + index" class="input-row">
          <el-input
            v-model="formData.scp_refs[index]"
            placeholder="/path/to/*.fa 或 /path/to/file.fasta"
            clearable
          >
            <template #append>
              <el-button @click="removeScpRef(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addScpRef">
          <el-icon><Plus /></el-icon> 添加 Reference 文件
        </el-button>
      </div>
    </div>

    <!-- Local Reference Files (upload mode) -->
    <div v-else-if="refSource === 'upload'">
      <div class="form-group">
        <label>或上传 Reference FASTA 文件</label>
        <div
          class="file-upload-zone"
          @click="triggerRefUpload"
        >
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 FASTA 文件</p>
        </div>
        <input
          ref="refUploadInput"
          type="file"
          accept=".fasta,.fa,.fna"
          style="display: none"
          @change="handleRefUpload"
        />
      </div>
    </div>

    <!-- NP Range -->
    <div class="form-group">
      <label>NP Range (可选)</label>
      <el-input
        v-model="formData.np_range"
        placeholder="例如: 0.9-1.0"
        clearable
      />
      <p class="hint">reads 的 nominal purity 范围</p>
    </div>

    <!-- RQ Range -->
    <div class="form-group">
      <label>RQ Range (可选)</label>
      <el-input
        v-model="formData.rq_range"
        placeholder="例如: 0.8-1.0"
        clearable
      />
      <p class="hint">reads 的质量范围</p>
    </div>

    <!-- Short Alignment -->
    <div class="form-group">
      <label>Short Alignment (可选)</label>
      <el-radio-group v-model="formData.short_aln">
        <el-radio :label="0">不使用 Short Alignment</el-radio>
        <el-radio :label="1">使用 Short Alignment</el-radio>
      </el-radio-group>
      <p class="hint">
        Short Alignment 用于处理长度在 [30, 200] 范围内的Reference序列。<br>
        如果您的Reference序列长度 < 200，建议选择"使用 Short Alignment"<br>
        如果Reference序列长度 ≥ 200，通常不需要Short Alignment
      </p>
    </div>

    <!-- Force Regenerate -->
    <div class="form-group">
      <el-checkbox v-model="formData.force">强制重新生成指标文件</el-checkbox>
      <p class="hint">如果指标文件已存在，默认会跳过。勾选此选项将强制重新生成</p>
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
const uploadInput = ref<HTMLInputElement | null>(null)

// Reference FASTA source and upload
const refSource = ref<'local' | 'upload' | 'scp'>('local')
const scpRefConfig = ref({
  server: '',
  password: ''
})
const refUploadInput = ref<HTMLInputElement | null>(null)

const formData = reactive({
  bams: ['/path/to/file1.bam'],
  refs: ['/path/to/reference.fa'],
  scp_refs: [] as string[],
  np_range: '',
  rq_range: '',
  short_aln: 0,
  force: false
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

const addRef = () => {
  formData.refs.push('/path/to/reference.fa')
}

const removeRef = (index: number) => {
  if (formData.refs.length > 1) {
    formData.refs.splice(index, 1)
  }
}

// SCP Reference file helpers
const addScpRef = () => {
  formData.scp_refs.push('/path/to/*.fasta 或 /path/to/file.fa')
}

const removeScpRef = (index: number) => {
  if (formData.scp_refs.length > 1) {
    formData.scp_refs.splice(index, 1)
  }
}

// Reference FASTA upload handlers
const triggerRefUpload = () => {
  refUploadInput.value?.click()
}

const handleRefUpload = async (event: Event) => {
  const target = event.target as HTMLInputElement
  const files = target.files
  if (!files) return

  for (let i = 0; i < files.length; i++) {
    const file = files[i]
    const validExtensions = ['.fasta', '.fa', '.fna']
    const hasValidExt = validExtensions.some(ext => file.name.toLowerCase().endsWith(ext))

    if (!hasValidExt) {
      ElMessage.warning(`只支持 .fasta, .fa, .fna 文件：${file.name}`)
      continue
    }

    try {
      const formDataObj = new FormData()
      formDataObj.append('file', file)

      const response = await axios.post('/api/files/upload', formDataObj)
      if (response.data.success) {
        formData.refs.push(response.data.data.local_path)
      }
    } catch (error: any) {
      ElMessage.error(`上传失败：${error.message || '未知错误'}`)
    }
  }

  if (refUploadInput.value) {
    refUploadInput.value.value = ''
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
      tool_name: 'reads-quality-hp',
      bams: formData.bams,
      refs: formData.refs,
      short_aln: formData.short_aln,
      force: formData.force
    }

    if (formData.np_range) {
      request.np_range = formData.np_range
    }
    if (formData.rq_range) {
      request.rq_range = formData.rq_range
    }

    // Add SSH config for remote BAM files
    if (fileSource.value === 'remote') {
      request.ssh_server = sshConfig.value.server
      request.ssh_password = sshConfig.value.password
    }

    // Add SCP reference files config
    if (refSource.value === 'scp') {
      if (!request.scp_refs) {
        request.scp_refs = []
      }
      request.scp_refs = formData.scp_refs

      if (scpRefConfig.value.server && scpRefConfig.value.password) {
        request.ssh_server = scpRefConfig.value.server
        request.ssh_password = scpRefConfig.value.password
      }
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
