<template>
  <div class="tool-form">
    <!-- BAM Basic Stat Form -->
    <BAMBasicStatForm v-if="selectedTool?.name === 'bam-basic-stat'" @execute="$emit('execute', $event)" />

    <!-- MSA View Form -->
    <MSAViewForm v-else-if="selectedTool?.name === 'msa-view'" @execute="$emit('execute', $event)" />

    <!-- Dynamic form for other tools -->
    <div v-else class="dynamic-form">
      <el-form :model="formData" label-width="120px" @submit.native.prevent="onSubmit">
        <template v-for="arg in argList" :key="arg.name">
          <el-form-item :label="arg.name" :prop="arg.name">
            <component :is="fieldComponent(arg)"
                       v-model="formData[arg.name]"
                       v-bind="fieldProps(arg)" />
          </el-form-item>
        </template>
        <el-form-item>
          <el-button type="primary" @click="onSubmit">执行工具</el-button>
        </el-form-item>
      </el-form>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import axios from 'axios'
import BAMBasicStatForm from './forms/BAMBasicStatForm.vue'
import MSAViewForm from './forms/MSAViewForm.vue'
import { ElMessage } from 'element-plus'

const props = defineProps<{ selectedTool: any }>()
const emit = defineEmits<{ execute: [data: any] }>()

const selectedTool = computed(() => props.selectedTool)

const argList = ref<Array<any>>([])
const formData = ref<Record<string, any>>({})

// Fetch argument schema when tool changes
watch(selectedTool, async (tool) => {
  if (!tool) return
  try {
    const resp = await axios.get(`/api/tools/${tool.name}/info`)
    if (resp.data.success) {
      argList.value = resp.data.data.arguments || []
      // Initialize formData with defaults
      const init: Record<string, any> = {}
      argList.value.forEach(arg => {
        if (arg.default !== undefined) init[arg.name] = arg.default
        else if (arg.type === 'boolean') init[arg.name] = false
        else init[arg.name] = ''
      })
      formData.value = init
    }
  } catch (e: any) {
    ElMessage.error('获取工具信息失败')
  }
})

const fieldComponent = (arg: any) => {
  switch (arg.type) {
    case 'select':
      return 'el-select'
    case 'number':
      return 'el-input-number'
    case 'file':
    case 'file_output':
      return 'el-input'
    default:
      return 'el-input'
  }
}

const fieldProps = (arg: any) => {
  const props: Record<string, any> = {}
  if (arg.type === 'select') {
    props.placeholder = arg.placeholder || ''
    props.filterable = true
  }
  if (arg.type === 'file' || arg.type === 'file_output') {
    props.placeholder = arg.placeholder || ''
    props.type = 'text'
  }
  if (arg.type === 'number') {
    props.placeholder = arg.placeholder || ''
    props.min = arg.min || undefined
    props.max = arg.max || undefined
  }
  if (arg.options) {
    props.options = arg.options
  }
  return props
}

const onSubmit = () => {
  emit('execute', { tool_name: selectedTool.value?.name, args: formData.value })
}
</script>

<style scoped>
.dynamic-form {
  background: #fff;
  padding: 30px;
  border-radius: 8px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.1);
}
</style>
