import axios from 'axios'

// API client configuration
const api = axios.create({
  baseURL: '/api',
  timeout: 60000,
  headers: {
    'Content-Type': 'application/json'
  }
})

// Tool types
export interface Tool {
  name: string
  module_path: string
}

export interface ToolInfo {
  name: string
  module_path: string
  description?: string
  arguments: any[]
}

export interface ToolExecutionResult {
  success: boolean
  tool_name: string
  stdout: string
  stderr: string
  exit_code: number
}

// API functions
export const toolsApi = {
  /** Get list of all available tools */
  async getToolsList(): Promise<Tool[]> {
    const response = await api.get('/tools/list')
    return response.data.success ? response.data.data.tools : []
  },

  /** Get information about a specific tool */
  async getToolInfo(toolName: string): Promise<ToolInfo> {
    const response = await api.get(`/tools/${toolName}/info`)
    return response.data.success ? response.data.data : null as any
  },

  /** Execute a CLI tool */
  async executeTool(toolName: string, args: any): Promise<ToolExecutionResult> {
    const response = await api.post(
      `/tools/${toolName}/execute`,
      { tool_name: toolName, args }
    )
    return response.data.success ? response.data.data : null as any
  }
}

export default api
