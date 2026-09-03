import axios from 'axios'

// API client configuration
const api = axios.create({
  baseURL: '/api',
  timeout: 60000,
  headers: {
    'Content-Type': 'application/json'
  }
})

// Add request interceptor for debugging
api.interceptors.request.use(
  (config) => {
    console.log('[Axios Request]', {
      method: config.method?.toUpperCase(),
      url: config.url,
      baseURL: config.baseURL,
      fullURL: config.url ? `${config.baseURL}${config.url}` : null,
      params: config.params,
      data: config.data,
      headers: config.headers
    })
    return config
  },
  (error) => {
    console.error('[Axios Request Error]', error)
    return Promise.reject(error)
  }
)

// Add response interceptor
api.interceptors.response.use(
  (response) => {
    console.log('[Axios Response]', {
      status: response.status,
      data: response.data
    })
    return response
  },
  (error) => {
    console.error('[Axios Response Error]', error)
    return Promise.reject(error)
  }
)

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
