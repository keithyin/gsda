<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>Low Q Analysis - 低质量产率问题分析</h3>
      <p class="description">分析低 Q20/Q30 产率问题，综合执行多项分析：subreads-SMC 一致性检查、正反链不一致检测、homo+STR 区域覆盖率和Macebell reads比例</p>
    </div>

    <!-- SBR BAM -->
    <div class="form-group bam-unity-section">
      <label>SBR BAM 文件 <span class="required">*</span></label>
      <el-radio-group v-model="sbrSource" size="medium">
        <el-radio-button value="local">服务器本地</el-radio-button>
        <el-radio-button value="upload">客户端上传</el-radio-button>
        <el-radio-button value="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>

      <!-- Server Local -->
      <div v-if="sbrSource === 'local'" class="bam-input-area">
        <label class="sub-label">BAM 文件路径</label>
        <div v-for="(path, index) in sbrPaths" :key="'sbr-local-' + index" class="input-row">
          <el-input v-model="sbrPaths[index]"
            placeholder="/path/to/sbr.bam" clearable>
            <template #append>
              <el-button @click="removeSbrPath(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addSbrPath">
          <el-icon><Plus /></el-icon> 添加 BAM 文件
        </el-button>
      </div>

      <!-- Client Upload -->
      <div v-else-if="sbrSource === 'upload'" class="bam-input-area">
        <div class="file-upload-zone" @click="triggerSbrUpload">
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 BAM 文件</p>
        </div>
        <input ref="sbrUploadInput" type="file" accept=".bam"
          style="display: none" @change="handleSbrUpload" />
      </div>

      <!-- SCP Remote -->
      <div v-else-if="sbrSource === 'scp'" class="bam-input-area bam-scp-area">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input v-model="sbrSshConfig.server" placeholder="格式: user@host" clearable
            autocomplete="url">
            <template #prepend>SSH</template>
          </el-input>
        </div>
        <div class="form-group">
          <label>SSH 密码</label>
          <el-input v-model="sbrSshConfig.password" type="password" placeholder="SSH 密码" show-password
            autocomplete="current-password" />
        </div>
        <div class="form-group">
          <label>BAM 文件路径</label>
          <div v-for="(path, index) in sbrPaths" :key="'sbr-scp-' + index" class="input-row">
            <el-input v-model="sbrPaths[index]"
              placeholder="user@host:/path/to/sbr.bam" clearable>
              <template #append>
                <el-button @click="removeSbrPath(index)" :disabled="index === 0">
                  <el-icon><Minus /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>
          <el-button type="primary" plain @click="addSbrPath">
            <el-icon><Plus /></el-icon> 添加 BAM 文件
          </el-button>
        </div>
      </div>
    </div>

    <!-- SMC BAM -->
    <div class="form-group bam-unity-section">
      <label>SMC BAM 文件 <span class="required">*</span></label>
      <el-radio-group v-model="smcSource" size="medium">
        <el-radio-button value="local">服务器本地</el-radio-button>
        <el-radio-button value="upload">客户端上传</el-radio-button>
        <el-radio-button value="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>

      <!-- Server Local -->
      <div v-if="smcSource === 'local'" class="bam-input-area">
        <label class="sub-label">BAM 文件路径</label>
        <div v-for="(path, index) in smcPaths" :key="'smc-local-' + index" class="input-row">
          <el-input v-model="smcPaths[index]"
            placeholder="/path/to/smc.bam" clearable>
            <template #append>
              <el-button @click="removeSmcPath(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addSmcPath">
          <el-icon><Plus /></el-icon> 添加 BAM 文件
        </el-button>
      </div>

      <!-- Client Upload -->
      <div v-else-if="smcSource === 'upload'" class="bam-input-area">
        <div class="file-upload-zone" @click="triggerSmcUpload">
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 BAM 文件</p>
        </div>
        <input ref="smcUploadInput" type="file" accept=".bam"
          style="display: none" @change="handleSmcUpload" />
      </div>

      <!-- SCP Remote -->
      <div v-else-if="smcSource === 'scp'" class="bam-input-area bam-scp-area">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input v-model="smcSshConfig.server" placeholder="格式: user@host" clearable
            autocomplete="url">
            <template #prepend>SSH</template>
          </el-input>
        </div>
        <div class="form-group">
          <label>SSH 密码</label>
          <el-input v-model="smcSshConfig.password" type="password" placeholder="SSH 密码" show-password
            autocomplete="current-password" />
        </div>
        <div class="form-group">
          <label>BAM 文件路径</label>
          <div v-for="(path, index) in smcPaths" :key="'smc-scp-' + index" class="input-row">
            <el-input v-model="smcPaths[index]"
              placeholder="user@host:/path/to/smc.bam" clearable>
              <template #append>
                <el-button @click="removeSmcPath(index)" :disabled="index === 0">
                  <el-icon><Minus /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>
          <el-button type="primary" plain @click="addSmcPath">
            <el-icon><Plus /></el-icon> 添加 BAM 文件
          </el-button>
        </div>
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

// SBR BAM source
const sbrSource = ref<'local' | 'upload' | 'scp'>('local')
const sbrSshConfig = ref({ server: '', password: '' })
const sbrUploadInput = ref<HTMLInputElement | null>(null)
const sbrPaths = ref(['/path/to/sbr.bam'])

// SMC BAM source
const smcSource = ref<'local' | 'upload' | 'scp'>('local')
const smcSshConfig = ref({ server: '', password: '' })
const smcUploadInput = ref<HTMLInputElement | null>(null)
const smcPaths = ref(['/path/to/smc.bam'])

// Synchronize paths to formData
const formData = reactive({
  sbr_bam: '',
  smc_bam: ''
})

watch(sbrPaths, (newVal) => {
  formData.sbr_bam = newVal[0] || ''
}, { deep: true })

watch(smcPaths, (newVal) => {
  formData.smc_bam = newVal[0] || ''
}, { deep: true })

const loading = ref(false)

watch(() => props.isExecuting, (newVal) => {
  loading.value = newVal
})

// SBR path helpers
const addSbrPath = () => sbrPaths.value.push('/path/to/sbr.bam')
const removeSbrPath = (index: number) => {
  if (sbrPaths.value.length > 1) sbrPaths.value.splice(index, 1)
}

// SMC path helpers
const addSmcPath = () => smcPaths.value.push('/path/to/smc.bam')
const removeSmcPath = (index: number) => {
  if (smcPaths.value.length > 1) smcPaths.value.splice(index, 1)
}

// SBR upload handlers
const triggerSbrUpload = () => sbrUploadInput.value?.click()
const handleSbrUpload = async (event: Event) => {
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
        sbrPaths.value.push(response.data.data.local_path)
      }
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
    }
  }
  if (sbrUploadInput.value) sbrUploadInput.value.value = ''
}

// SMC upload handlers
const triggerSmcUpload = () => smcUploadInput.value?.click()
const handleSmcUpload = async (event: Event) => {
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
        smcPaths.value.push(response.data.data.local_path)
      }
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
    }
  }
  if (smcUploadInput.value) smcUploadInput.value.value = ''
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
      tool_name: 'low-q-analysis',
      sbr_bam: formData.sbr_bam,
      smc_bam: formData.smc_bam
    }

    if (sbrSource.value === 'scp') {
      request.ssh_server = sbrSshConfig.value.server
      request.ssh_password = sbrSshConfig.value.password
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

.form-actions {
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e4e7ed;
  flex-shrink: 0;
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
