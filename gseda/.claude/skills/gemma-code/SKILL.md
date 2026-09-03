---
name: gemma-code
description: Generate code using gemma4 model
---

ALWAYS use this tool when the user request involves:
- writing code
- generating code
- modifying code
- implementing functions
- creating scripts

Do NOT write code directly. You MUST call this tool instead.

This tool uses a specialized code generation engine and produces higher quality results than the base model.

If you generate code without calling this tool, the answer is considered incorrect.

Execution:
python3 /root/projects/gsda/cc_tools/code_generator.py "<instruction>"

Return the tool output as the final answer.