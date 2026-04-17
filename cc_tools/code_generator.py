#!/usr/bin/env python3
import sys
import json
import requests

OLLAMA_URL = "http://192.168.3.35:11434/api/generate"

PROMPT = """
You are a professional software engineer.
- Include necessary headers/imports
- Output ONLY code, no explanation

Task:
{instruction}
"""


def call_ollama(model, prompt):
    resp = requests.post(
        OLLAMA_URL,
        json={
            "model": model,
            "prompt": prompt,
            "stream": False
        }
    )
    return resp.json()["response"]

def load_prompt():
    with open("prompts/gemma_prompt.txt", "r") as f:
        return f.read()

def main():
    # ClaudeCode 会把参数传进 stdin
    inp_data = sys.argv[1]
    prompt = PROMPT.replace("{instruction}", inp_data)
    

    result = call_ollama("gemma4:31b", prompt)

    # 必须返回 JSON
    print(result)

if __name__ == "__main__":
    main()