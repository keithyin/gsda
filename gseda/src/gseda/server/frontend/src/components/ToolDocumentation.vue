<template>
  <div v-if="doc" class="documentation-content">
    <el-tabs v-model="activeTab">
      <!-- 功能概述 -->
      <el-tab-pane label="功能概述" name="overview">
        <el-card shadow="never" class="overview-card">
          <template #header>
            <span class="card-title">{{ doc.title }}</span>
          </template>
          <p class="overview-text">{{ doc.overview }}</p>
        </el-card>
      </el-tab-pane>

      <!-- 输入参数详解 -->
      <el-tab-pane label="输入参数" name="parameters">
        <el-table :data="doc.parameters" border size="default" class="params-table">
          <el-table-column prop="name" label="参数名" width="140" />
          <el-table-column prop="cli_flag" label="前端字段" width="160">
            <template #default="{ row }">
              <code v-if="row.cli_flag === 'positional'">位置参数</code>
              <code v-else>{{ row.cli_flag }}</code>
            </template>
          </el-table-column>
          <el-table-column prop="type" label="类型" width="80" align="center">
            <template #default="{ row }">
              <el-tag size="small" type="info">{{ row.type }}</el-tag>
            </template>
          </el-table-column>
          <el-table-column prop="required" label="必需" width="70" align="center">
            <template #default="{ row }">
              <el-tag :type="row.required ? 'danger' : 'info'" size="small">
                {{ row.required ? '是' : '否' }}
              </el-tag>
            </template>
          </el-table-column>
          <el-table-column label="默认值" width="90" align="center">
            <template #default="{ row }">
              <code>{{ row.default ?? '-' }}</code>
            </template>
          </el-table-column>
          <el-table-column prop="description" label="描述" min-width="300" />
        </el-table>
      </el-tab-pane>

      <!-- 使用示例 -->
      <el-tab-pane label="使用示例" name="examples">
        <div v-for="(ex, i) in doc.examples" :key="i" class="example-block">
          <h4 class="example-title">{{ ex.title }}</h4>
          <el-card shadow="never" class="steps-card">
            <ol class="steps-list">
              <li v-for="(step, j) in ex.steps" :key="j" class="step-item">{{ step }}</li>
            </ol>
          </el-card>
        </div>
      </el-tab-pane>

      <!-- 输出说明 -->
      <el-tab-pane label="输出说明" name="output">
        <p v-if="doc.output?.description" class="output-desc">{{ doc.output.description }}</p>
        <div v-for="(section, i) in doc.output?.sections" :key="i" class="output-section">
          <h4 class="output-title">{{ section.name }}</h4>
          <p class="output-text">{{ section.content }}</p>
        </div>
      </el-tab-pane>
    </el-tabs>
  </div>
  <el-empty v-else description="文档加载中..." :image-size="60" />
</template>

<script setup lang="ts">
import { ref } from 'vue'

defineProps<{
  doc: any
}>()

const activeTab = ref('overview')
</script>

<style scoped>
.documentation-content {
  max-height: 70vh;
  overflow-y: auto;
}

.overview-card {
  margin-bottom: 0;
}

.card-title {
  font-size: 16px;
  font-weight: 600;
  color: #303133;
}

.overview-text {
  margin: 0;
  line-height: 1.8;
  color: #606266;
  font-size: 14px;
}

.params-table {
  margin-top: 10px;
}

.params-table :deep(.el-table th) {
  background: #f5f7fa;
  color: #606266;
  font-weight: 600;
}

.params-table code {
  background: #f5f7fa;
  padding: 2px 6px;
  border-radius: 3px;
  font-size: 13px;
}

.example-block {
  margin-bottom: 20px;
}

.example-block:last-child {
  margin-bottom: 0;
}

.example-title {
  margin: 0 0 10px 0;
  font-size: 15px;
  color: #303133;
  font-weight: 600;
}

.steps-card :deep(.el-card__body) {
  padding: 12px 16px;
}

.steps-list {
  margin: 0;
  padding-left: 20px;
  line-height: 2;
  color: #606266;
  font-size: 14px;
}

.steps-list .step-item {
  margin-bottom: 8px;
  white-space: pre-wrap;
}

.output-desc {
  margin: 0 0 16px 0;
  color: #606266;
  font-size: 14px;
}

.output-section {
  margin-bottom: 16px;
  padding: 12px 16px;
  background: #fafbfc;
  border-radius: 4px;
  border-left: 3px solid #409eff;
}

.output-section:last-child {
  margin-bottom: 0;
}

.output-title {
  margin: 0 0 6px 0;
  font-size: 14px;
  color: #303133;
  font-weight: 600;
}

.output-text {
  margin: 0;
  color: #606266;
  font-size: 13px;
  line-height: 1.6;
}

:deep(.el-tabs__header) {
  margin: 0 0 16px 0;
}

:deep(.el-tabs__content) {
  padding: 0;
}
</style>
