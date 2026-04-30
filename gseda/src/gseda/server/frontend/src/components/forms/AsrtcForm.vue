<template>
  <div class="form-container" :class="{ 'is-loading': isExecuting }">
    <div class="form-section">
      <h3>ASRTC - ASmRTP Consistency</h3>
      <p class="description">比较 subreads 和 SMC 的一致性分析，输出 asrtc.txt 文件</p>
    </div>

    <!-- ref_fa -->
    <div class="form-group bam-unity-section">
      <label>参考序列 FASTA <span class="required">*</span></label>
      <el-radio-group v-model="refFaSource" size="medium">
        <el-radio-button value="local">服务器本地</el-radio-button>
        <el-radio-button value="upload">客户端上传</el-radio-button>
        <el-radio-button value="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>

      <div v-if="refFaSource === 'local'" class="bam-input-area">
        <label class="sub-label">文件路径</label>
        <div v-for="(path, index) in refFaPaths" :key="'refFa-local-' + index" class="input-row">
          <el-input v-model="refFaPaths[index]" placeholder="/path/to/ref.fa" clearable>
            <template #append>
              <el-button @click="removeRefFaPath(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addRefFaPath">
          <el-icon><Plus /></el-icon> 添加文件
        </el-button>
      </div>

      <div v-else-if="refFaSource === 'upload'" class="bam-input-area">
        <div class="file-upload-zone" @click="triggerRefFaUpload">
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 FASTA 文件</p>
        </div>
        <input ref="refFaUploadInput" type="file" accept=".fasta,.fa"
          style="display: none" @change="handleRefFaUpload" />
      </div>

      <div v-else-if="refFaSource === 'scp'" class="bam-input-area bam-scp-area">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input v-model="refFaSshConfig.server" placeholder="格式: user@host" clearable autocomplete="url">
            <template #prepend>SSH</template>
          </el-input>
        </div>
        <div class="form-group">
          <label>SSH 密码</label>
          <el-input v-model="refFaSshConfig.password" type="password" placeholder="SSH 密码" show-password
            autocomplete="current-password" />
        </div>
        <div class="form-group">
          <label>文件路径</label>
          <div v-for="(path, index) in refFaPaths" :key="'refFa-scp-' + index" class="input-row">
            <el-input v-model="refFaPaths[index]" placeholder="user@host:/path/to/ref.fa" clearable>
              <template #append>
                <el-button @click="removeRefFaPath(index)" :disabled="index === 0">
                  <el-icon><Minus /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>
          <el-button type="primary" plain @click="addRefFaPath">
            <el-icon><Plus /></el-icon> 添加文件
          </el-button>
        </div>
      </div>
    </div>

    <!-- sbr -->
    <div class="form-group bam-unity-section">
      <label>Subreads BAM <span class="required">*</span></label>
      <el-radio-group v-model="sbrSource" size="medium">
        <el-radio-button value="local">服务器本地</el-radio-button>
        <el-radio-button value="upload">客户端上传</el-radio-button>
        <el-radio-button value="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>

      <div v-if="sbrSource === 'local'" class="bam-input-area">
        <label class="sub-label">文件路径</label>
        <div v-for="(path, index) in sbrPaths" :key="'sbr-local-' + index" class="input-row">
          <el-input v-model="sbrPaths[index]" placeholder="/path/to/subreads.bam" clearable>
            <template #append>
              <el-button @click="removeSbrPath(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addSbrPath">
          <el-icon><Plus /></el-icon> 添加文件
        </el-button>
      </div>

      <div v-else-if="sbrSource === 'upload'" class="bam-input-area">
        <div class="file-upload-zone" @click="triggerSbrUpload">
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 BAM 文件</p>
        </div>
        <input ref="sbrUploadInput" type="file" accept=".bam"
          style="display: none" @change="handleSbrUpload" />
      </div>

      <div v-else-if="sbrSource === 'scp'" class="bam-input-area bam-scp-area">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input v-model="sbrSshConfig.server" placeholder="格式: user@host" clearable autocomplete="url">
            <template #prepend>SSH</template>
          </el-input>
        </div>
        <div class="form-group">
          <label>SSH 密码</label>
          <el-input v-model="sbrSshConfig.password" type="password" placeholder="SSH 密码" show-password
            autocomplete="current-password" />
        </div>
        <div class="form-group">
          <label>文件路径</label>
          <div v-for="(path, index) in sbrPaths" :key="'sbr-scp-' + index" class="input-row">
            <el-input v-model="sbrPaths[index]" placeholder="user@host:/path/to/subreads.bam" clearable>
              <template #append>
                <el-button @click="removeSbrPath(index)" :disabled="index === 0">
                  <el-icon><Minus /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>
          <el-button type="primary" plain @click="addSbrPath">
            <el-icon><Plus /></el-icon> 添加文件
          </el-button>
        </div>
      </div>
    </div>

    <!-- smc -->
    <div class="form-group bam-unity-section">
      <label>SMC BAM/FASTA <span class="required">*</span></label>
      <el-radio-group v-model="smcSource" size="medium">
        <el-radio-button value="local">服务器本地</el-radio-button>
        <el-radio-button value="upload">客户端上传</el-radio-button>
        <el-radio-button value="scp">SCP 远程文件</el-radio-button>
      </el-radio-group>

      <div v-if="smcSource === 'local'" class="bam-input-area">
        <label class="sub-label">文件路径</label>
        <div v-for="(path, index) in smcPaths" :key="'smc-local-' + index" class="input-row">
          <el-input v-model="smcPaths[index]" placeholder="/path/to/smc.bam" clearable>
            <template #append>
              <el-button @click="removeSmcPath(index)" :disabled="index === 0">
                <el-icon><Minus /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <el-button type="primary" plain @click="addSmcPath">
          <el-icon><Plus /></el-icon> 添加文件
        </el-button>
      </div>

      <div v-else-if="smcSource === 'upload'" class="bam-input-area">
        <div class="file-upload-zone" @click="triggerSmcUpload">
          <el-icon :size="48" color="#909399"><Upload /></el-icon>
          <p style="margin-top: 10px; color: #909399;">点击上传 BAM/FASTA 文件</p>
        </div>
        <input ref="smcUploadInput" type="file" accept=".bam,.fasta,.fa,.fq,.fastq"
          style="display: none" @change="handleSmcUpload" />
      </div>

      <div v-else-if="smcSource === 'scp'" class="bam-input-area bam-scp-area">
        <div class="form-group">
          <label>SSH 服务器地址</label>
          <el-input v-model="smcSshConfig.server" placeholder="格式: user@host" clearable autocomplete="url">
            <template #prepend>SSH</template>
          </el-input>
        </div>
        <div class="form-group">
          <label>SSH 密码</label>
          <el-input v-model="smcSshConfig.password" type="password" placeholder="SSH 密码" show-password
            autocomplete="current-password" />
        </div>
        <div class="form-group">
          <label>文件路径</label>
          <div v-for="(path, index) in smcPaths" :key="'smc-scp-' + index" class="input-row">
            <el-input v-model="smcPaths[index]" placeholder="user@host:/path/to/smc.bam" clearable>
              <template #append>
                <el-button @click="removeSmcPath(index)" :disabled="index === 0">
                  <el-icon><Minus /></el-icon>
                </el-button>
              </template>
            </el-input>
          </div>
          <el-button type="primary" plain @click="addSmcPath">
            <el-icon><Plus /></el-icon> 添加文件
          </el-button>
        </div>
      </div>
    </div>

    <!-- prefix -->
    <div class="form-group">
      <label>输出前缀 <span class="required">*</span></label>
      <el-input v-model="formData.prefix" placeholder="输出前缀 (输出文件: {prefix}.asrtc.txt)" />
    </div>

    <!-- rq_range -->
    <div class="form-group">
      <label>RQ 范围过滤</label>
      <el-input v-model="formData.rq_range" placeholder="例如: 0.8:1.0" />
    </div>

    <!-- np_range -->
    <div class="form-group">
      <label>NP 范围过滤</label>
      <el-input v-model="formData.np_range" placeholder="例如: 1:3,5,7:9" />
    </div>

    <!-- Execute Button -->
    <div class="form-actions">
      <el-button type="primary" size="large" :loading="loading" @click="execute">
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

// ref_fa
const refFaSource = ref<'local' | 'upload' | 'scp'>('local')
const refFaSshConfig = ref({ server: '', password: '' })
const refFaUploadInput = ref<HTMLInputElement | null>(null)
const refFaPaths = ref(['/path/to/ref.fa'])

// sbr
const sbrSource = ref<'local' | 'upload' | 'scp'>('local')
const sbrSshConfig = ref({ server: '', password: '' })
const sbrUploadInput = ref<HTMLInputElement | null>(null)
const sbrPaths = ref(['/path/to/subreads.bam'])

// smc
const smcSource = ref<'local' | 'upload' | 'scp'>('local')
const smcSshConfig = ref({ server: '', password: '' })
const smcUploadInput = ref<HTMLInputElement | null>(null)
const smcPaths = ref(['/path/to/smc.bam'])

// formData
const formData = reactive({
  ref_fa: '',
  sbr: '',
  smc: '',
  prefix: '',
  rq_range: '',
  np_range: ''
})

watch(refFaPaths, (newVal) => { formData.ref_fa = newVal[0] || '' }, { deep: true })
watch(sbrPaths, (newVal) => { formData.sbr = newVal[0] || '' }, { deep: true })
watch(smcPaths, (newVal) => { formData.smc = newVal[0] || '' }, { deep: true })

const loading = ref(false)
watch(() => props.isExecuting, (newVal) => { loading.value = newVal })

// ref_fa helpers
const addRefFaPath = () => refFaPaths.value.push('/path/to/ref.fa')
const removeRefFaPath = (index: number) => { if (refFaPaths.value.length > 1) refFaPaths.value.splice(index, 1) }
const triggerRefFaUpload = () => refFaUploadInput.value?.click()
const handleRefFaUpload = async (event: Event) => {
  const target = event.target as HTMLInputElement
  const files = target.files
  if (!files) return
  for (const file of Array.from(files)) {
    const validExts = ['.fasta', '.fa']
    const ext = '.' + file.name.split('.').pop()?.toLowerCase()
    if (!validExts.includes(ext)) {
      ElMessage.warning(`只支持 .fasta, .fa 文件: ${file.name}`)
      continue
    }
    try {
      const fd = new FormData()
      fd.append('file', file)
      const response = await axios.post('/api/files/upload', fd)
      if (response.data.success) refFaPaths.value.push(response.data.data.local_path)
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
    }
  }
  if (refFaUploadInput.value) refFaUploadInput.value.value = ''
}

// sbr helpers
const addSbrPath = () => sbrPaths.value.push('/path/to/subreads.bam')
const removeSbrPath = (index: number) => { if (sbrPaths.value.length > 1) sbrPaths.value.splice(index, 1) }
const triggerSbrUpload = () => sbrUploadInput.value?.click()
const handleSbrUpload = async (event: Event) => {
  const target = event.target as HTMLInputElement
  const files = target.files
  if (!files) return
  for (const file of Array.from(files)) {
    if (!file.name.toLowerCase().endsWith('.bam')) {
      ElMessage.warning(`只支持 .bam 文件: ${file.name}`)
      continue
    }
    try {
      const fd = new FormData()
      fd.append('file', file)
      const response = await axios.post('/api/files/upload', fd)
      if (response.data.success) sbrPaths.value.push(response.data.data.local_path)
    } catch (error: any) {
      ElMessage.error(`上传失败: ${error.message || '未知错误'}`)
    }
  }
  if (sbrUploadInput.value) sbrUploadInput.value.value = ''
}

// smc helpers
const addSmcPath = () => smcPaths.value.push('/path/to/smc.bam')
const removeSmcPath = (index: number) => { if (smcPaths.value.length > 1) smcPaths.value.splice(index, 1) }
const triggerSmcUpload = () => smcUploadInput.value?.click()
const handleSmcUpload = async (event: Event) => {
  const target = event.target as HTMLInputElement
  const files = target.files
  if (!files) return
  for (const file of Array.from(files)) {
    const ext = '.' + file.name.split('.').pop()?.toLowerCase()
    const validExts = ['.bam', '.fasta', '.fa', '.fq', '.fastq']
    if (!validExts.includes(ext)) {
      ElMessage.warning(`只支持 .bam, .fasta, .fa, .fq, .fastq 文件: ${file.name}`)
      continue
    }
    try {
      const fd = new FormData()
      fd.append('file', file)
      const response = await axios.post('/api/files/upload', fd)
      if (response.data.success) smcPaths.value.push(response.data.data.local_path)
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
    if (actions) actions.scrollIntoView({ behavior: 'smooth', block: 'center' })
  }
}

const execute = async () => {
  await nextTick()
  scrollToActions()
  loading.value = true
  try {
    const request: any = {
      tool_name: 'asrtc',
      ref_fa: formData.ref_fa,
      sbr: formData.sbr,
      smc: formData.smc,
      prefix: formData.prefix,
    }
    if (formData.rq_range) request.rq_range = formData.rq_range
    if (formData.np_range) request.np_range = formData.np_range

    // Carry over SSH config from whichever source selected SCP
    if (refFaSource.value === 'scp') {
      request.ssh_server = refFaSshConfig.value.server
      request.ssh_password = refFaSshConfig.value.password
    }
    if (sbrSource.value === 'scp') {
      request.ssh_server = sbrSshConfig.value.server
      request.ssh_password = sbrSshConfig.value.password
    }
    if (smcSource.value === 'scp') {
      request.ssh_server = smcSshConfig.value.server
      request.ssh_password = smcSshConfig.value.password
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
  top: 0; left: 0; right: 0; bottom: 0;
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
