// TypeScript type definitions for GSEDA Frontend

export interface Tool {
  name: string
  module_path: string
}

export interface ToolInfo {
  name: string
  module_path: string
  description?: string
  arguments: ArgumentSchema[]
}

export interface ArgumentSchema {
  name: string
  type: 'text' | 'number' | 'select' | 'file' | 'textarea' | 'file_output'
  required?: boolean
  multiple?: boolean
  placeholder?: string
  default?: any
  options?: string[]
}

export interface ToolExecutionRequest {
  tool_name: string
  args: Record<string, any>
}

export interface ToolExecutionResponse {
  success: boolean
  tool_name: string
  stdout: string
  stderr: string
  exit_code: number
  result_data?: Record<string, any>
}

export interface APIResponse<T = any> {
  success: boolean
  data?: T
  error?: string
  message?: string
}
