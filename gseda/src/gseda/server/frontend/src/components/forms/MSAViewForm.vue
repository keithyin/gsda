<template>
  <div class="form-container">
    <div class="form-section">
      <h3>MSA View - 多重序列比对视图</h3>
      <p class="description">将 BAM 文件中的比对结果转换为 MSA（多重序列比对）格式，支持生成 FASTA 文件和可视化图片。</p>
    </div>

    <!-- BAM File -->
    <div class="form-group">
      <label>Alignment BAM</label>
      <el-input
        v-model="formData.bam"
        placeholder="/path/to/aligned.bam"
        clearable
      />
    </div>

    <!-- Reference FASTA -->
    <div class="form-group">
      <label>Reference FASTA</label>
      <el-input
        v-model="formData.ref_fasta"
        placeholder="/path/to/reference.fa"
        clearable
      >
        <template #prefix>
        </template>
      </el-input>
    </div>

    <!-- Reference Name -->
    <div class="form-group">
      <label>Reference Name</label>
      <el-input
        v-model="formData.ref_name"
        placeholder="例如：contig_001"
        clearable
      />
    </div>

    <!-- Reference Start (optional) -->
    <div class="form-group">
      <label>Reference Start (可选)</label>
      <el-input
        v-model.number="formData.start"
        type="number"
        placeholder="起始位置（留空则自动检测）"
      />
    </div>

    <!-- Reference End (optional) -->
    <div class="form-group">
      <label>Reference End (可选)</label>
      <el-input
        v-model.number="formData.end"
        type="number"
        placeholder="结束位置（留空则自动检测）"
      />
    </div>

    <!-- Output FASTA -->
    <div class="form-group">
      <label>输出 FASTA 文件 (可选)</label>
      <el-input
        v-model="formData.o_fasta"
        placeholder="/path/to/output.fasta"
        clearable
      >
        <template #prefix>
        </template>
      </el-input>
    </div>

    <!-- Output Picture -->
    <div class="form-group">
      <label>输出图片 (可选)</label>
      <el-input
        v-model="formData.o_pic"
        placeholder="/path/to/output.png"
        clearable
      >
        <template #prefix>
        </template>
      </el-input>
    </div>

    <!-- Execute Button -->
    <div class="form-actions">
      <el-button type="primary" size="large" :loading="loading" @click="execute">
        {{ loading ? '生成中...' : '生成 MSA' }}
      </el-button>
    </div>
  </div>
</template>

<script setup lang="ts">
import { ref, reactive, onMounted } from 'vue'

const emit = defineEmits<{
  execute: [data: any]
}>()

const formData = reactive({
  bam: '',
  ref_fasta: '',
  ref_name: '',
  start: null as number | null,
  end: null as number | null,
  o_fasta: null as string | null,
  o_pic: null as string | null
})

const loading = ref(false)

onMounted(() => {
  // Watch for changes and update parent
  watch(formData, (newVal) => {
    emit('execute', { tool_name: 'msa-view', args: newVal })
  }, { deep: true })
})

const execute = async () => {
  emit('execute', { tool_name: 'msa-view', args: { ...formData } })
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
  margin-bottom: 20px;
}

.form-group label {
  display: block;
  margin-bottom: 8px;
  color: #606266;
  font-weight: 500;
}

.input-with-icon {
  position: relative;
}

.el-input__prefix {
  padding-left: 8px;
  color: #c0c4cc;
}

.form-actions {
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e4e7ed;
}
</style>
