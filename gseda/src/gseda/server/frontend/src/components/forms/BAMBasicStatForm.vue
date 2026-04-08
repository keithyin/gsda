<template>
  <div class="form-container">
    <div class="form-section">
      <h3>BAM Basic Stat - BAM 文件基本统计</h3>
      <p class="description">统计 BAM 文件中 reads 的基本质量信息，包括 channel 分布、测序长度分布、质量分数统计等。</p>
    </div>

    <!-- BAM Files -->
    <div class="form-group">
      <label>
        BAM 文件路径
        <span class="required">*</span>
      </label>
      <div v-for="(path, index) in formData.bams" :key="index" class="input-row">
        <el-input
          v-model="formData.bams[index]"
          placeholder="/path/to/file.bam"
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
      <el-button type="primary" size="large" :loading="loading" @click="execute">
        {{ loading ? '分析中...' : '开始分析' }}
      </el-button>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, reactive, watch } from 'vue'
import { Plus, Minus } from '@element-plus/icons-vue'

const emit = defineEmits<{
  execute: [data: any]
}>()

const formData = reactive({
  bams: ['/path/to/file1.bam'],
  channel_tag: 'ch',
  min_rq: null as number | null
})

const loading = ref(false)

// Watch for changes and update parent
watch(formData, (newVal) => {
  emit('execute', { tool_name: 'bam-basic-stat', args: newVal })
}, { deep: true })

const addBam = () => {
  formData.bams.push('/path/to/file.bam')
}

const removeBam = (index: number) => {
  if (formData.bams.length > 1) {
    formData.bams.splice(index, 1)
  }
}

const execute = async () => {
  emit('execute', { tool_name: 'bam-basic-stat', args: { ...formData } })
}
</script>

<style scoped>
.form-container {
  background: #fff;
  padding: 25px;
  border-radius: 8px;
  box-shadow: 0 2px 12px rgba(0, 0, 0, 0.1);
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

.form-actions {
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e4e7ed;
}
</style>
