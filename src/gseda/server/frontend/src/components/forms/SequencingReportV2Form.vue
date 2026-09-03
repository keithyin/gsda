<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>Sequencing Report V2 - 综合测序质量报告</h3>
      <p class="description">基于 BAM 和 Reference 生成综合测序质量报告（含比对统计和基础统计），或合并已有的 Fact/Basic CSV 指标文件</p>
    </div>

    <!-- Mode Toggle -->
    <div class="form-group">
      <label>运行模式</label>
      <el-radio-group v-model="mode" size="medium">
        <el-radio-button value="normal">正常分析</el-radio-button>
        <el-radio-button value="merge">合并已有CSV</el-radio-button>
      </el-radio-group>
    </div>

    <!-- Normal Mode: BAM Source Toggle -->
    <div v-if="mode === 'normal'">
      <div class="form-group">
        <label>BAM 文件来源</label>
        <el-radio-group v-model="fileSource" size="medium">
          <el-radio-button value="local">本地文件</el-radio-button>
          <el-radio-button value="remote">远程文件 (SCP)</el-radio-button>
        </el-radio-group>
      </div>

      <!-- Remote BAM File Inputs -->
      <div v-if="fileSource === 'remote'" class="remote-config">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input
            v-model="sshConfig.server"
            placeholder="格式: user@host"
            clearable
            autocomplete="url"
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
            autocomplete="current-password"
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

      <!-- Local BAM File Inputs -->
      <div v-else>
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

        <div class="form-group">
          <label>或上传 BAM 文件</label>
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
    </div>

    <!-- Merge Mode: CSV file inputs -->
    <div v-if="mode === 'merge'">
      <div class="form-group">
        <label>Fact Metric CSV 文件 (可选)</label>
        <div v-for="(path, index) in formData.fact_csvs" :key="'fact-' + index" class="input-row">
          <el-input
            v-model="formData.fact_csvs[index]"
            placeholder="/path/to/*.csv"
            clearable
          >
            <template #append>
              <el-button @click="removeFactCsv(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addFactCsv">
          <el-icon><Plus /></el-icon> 添加 Fact CSV 文件
        </el-button>
      </div>

      <div class="form-group">
        <label>Basic CSV 文件 (可选)</label>
        <div v-for="(path, index) in formData.basic_csvs" :key="'basic-' + index" class="input-row">
          <el-input
            v-model="formData.basic_csvs[index]"
            placeholder="/path/to/*.csv"
            clearable
          >
            <template #append>
              <el-button @click="removeBasicCsv(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addBasicCsv">
          <el-icon><Plus /></el-icon> 添加 Basic CSV 文件
        </el-button>
      </div>
    </div>

    <!-- Reference FASTA Files (normal mode only) -->
    <div v-if="mode === 'normal'">
      <div class="form-group">
        <label>Reference FASTA 文件 (可选)</label>
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
        <p class="hint">如果不提供，则不输出对齐相关的指标</p>
      </div>

      <!-- Reference FASTA Source Toggle -->
      <div class="form-group">
        <label>Reference FASTA 文件来源</label>
        <el-radio-group v-model="refSource" size="medium">
          <el-radio-button value="local">服务器本地</el-radio-button>
          <el-radio-button value="upload">客户端上传</el-radio-button>
          <el-radio-button value="scp">SCP 远程文件</el-radio-button>
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
            autocomplete="url"
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
            autocomplete="current-password"
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
    </div>

    <!-- Short Alignment -->
    <div class="form-group">
      <label>Short Alignment (可选)</label>
      <el-radio-group v-model="formData.short_aln">
        <el-radio :value="0">不使用 Short Alignment</el-radio>
        <el-radio :value="1">使用 Short Alignment</el-radio>
      </el-radio-group>
      <p class="hint">
        Short Alignment 用于处理长度在 [30, 200] 范围内的Reference序列。<br>
        如果您的Reference序列长度 < 200，建议选择"使用 Short Alignment"<br>
        如果Reference序列长度 ≥ 200，通常不需要Short Alignment
      </p>
    </div>

    <!-- Disable Basic Stat -->
    <div class="form-group">
      <el-checkbox v-model="formData.disable_basic_stat">跳过基础统计计算</el-checkbox>
      <p class="hint">跳过 bam_basic 统计，仅输出对齐相关指标</p>
    </div>

    <!-- Disable Align Stat -->
    <div class="form-group">
      <el-checkbox v-model="formData.disable_align_stat">跳过对齐统计计算</el-checkbox>
      <p class="hint">跳过 gsmm2-aligned-metric 对齐统计，仅输出基础统计</p>
    </div>

    <!-- np_range -->
    <div class="form-group">
      <label>np 范围过滤 (可选)</label>
      <el-input
        v-model="formData.np_range"
        placeholder="例如: 1:3,5,7:9 表示 [[1,3], [5,5], [7,9]]"
        clearable
      />
      <p class="hint">仅对包含 np 字段的 BAM 文件有效</p>
    </div>

    <!-- rq_range -->
    <div class="form-group">
      <label>rq 范围过滤 (可选)</label>
      <el-input
        v-model="formData.rq_range"
        placeholder="例如: 0.9:1.1 表示 0.9≤rq≤1.1"
        clearable
      />
      <p class="hint">仅对包含 rq 字段的 BAM 文件有效</p>
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

const mode = ref<'normal' | 'merge'>('normal')
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
  bams: ['/path/to/file1.bam'] as string[],
  refs: [] as string[],
  scp_refs: [] as string[],
  fact_csvs: ['/path/to/fact1.csv'] as string[],
  basic_csvs: ['/path/to/basic1.csv'] as string[],
  short_aln: 0 as number,
  disable_basic_stat: false as boolean,
  disable_align_stat: false as boolean,
  np_range: '' as string,
  rq_range: '' as string,
  force: true as boolean
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
  if (formData.refs.length > 0) {
    formData.refs.splice(index, 1)
  }
}

const addFactCsv = () => {
  formData.fact_csvs.push('/path/to/fact.csv')
}

const removeFactCsv = (index: number) => {
  if (formData.fact_csvs.length > 1) {
    formData.fact_csvs.splice(index, 1)
  }
}

const addBasicCsv = () => {
  formData.basic_csvs.push('/path/to/basic.csv')
}

const removeBasicCsv = (index: number) => {
  if (formData.basic_csvs.length > 1) {
    formData.basic_csvs.splice(index, 1)
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
      tool_name: 'sequencing-report-v2'
    }

    if (mode.value === 'normal') {
      request.bams = formData.bams
      request.refs = formData.refs.length > 0 ? formData.refs : undefined
    } else {
      request.fact_csvs = formData.fact_csvs.length > 0 ? formData.fact_csvs : undefined
      request.basic_csvs = formData.basic_csvs.length > 0 ? formData.basic_csvs : undefined
    }

    if (formData.short_aln !== 0) {
      request.short_aln = formData.short_aln
    }
    if (formData.disable_basic_stat) {
      request.disable_basic_stat = formData.disable_basic_stat
    }
    if (formData.disable_align_stat) {
      request.disable_align_stat = formData.disable_align_stat
    }
    if (formData.np_range) {
      request.np_range = formData.np_range
    }
    if (formData.rq_range) {
      request.rq_range = formData.rq_range
    }
    if (formData.force) {
      request.force = formData.force
    }

    // Add SSH config for remote BAM files
    if (mode.value === 'normal' && fileSource.value === 'remote') {
      request.ssh_server = sshConfig.value.server
      request.ssh_password = sshConfig.value.password
    }

    // Add SCP reference files config
    if (mode.value === 'normal' && refSource.value === 'scp') {
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
