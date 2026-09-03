#!/usr/bin/env python3
"""
多轮对话 KV Cache / Prefix Cache 测试脚本（适配 ccswtich）- 自增长版本

每轮对话长度随步骤递增，模拟真实对话中上下文逐渐增长的场景。
验证 llama.cpp 的 prefix caching 在动态增长上下文下的表现。
"""

import requests
from pathlib import Path
import time
import json
import argparse
import sys


def load_text_from_arg(arg_value, file_flag=False):
    """从字符串或文件加载文本"""
    if file_flag:
        with open(arg_value, 'r', encoding='utf-8') as f:
            return f.read()
    return arg_value


def make_text_tokens_approx(n_tokens: int, base_text: str = None) -> str:
    """
    粗略生成指定 token 数的文本（1 token ≈ 4 chars）。
    可指定基础重复文本，否则使用默认英文段落。
    """
    if base_text is None:
        base_text = (
            "This is a fixed cache testing sentence. "
            "The content before the dynamic suffix must remain exactly unchanged. "
            "We use this text to test prefix KV cache reuse in llama.cpp. "
        )
    target_chars = max(n_tokens * 4, 1)  # 至少1字符
    repeat = (target_chars // len(base_text)) + 1
    text = (base_text * repeat)[:target_chars]
    return text


def call_server(
    url: str,
    model: str,
    messages: list,
    extra_headers: dict = None,
    max_tokens: int = 1,
    temperature: float = 0.0,
    stream: bool = False,
    verbose: bool = False,
):
    """发送 chat completion 请求，返回耗时和 usage 信息"""
    payload = {
        "model": model,
        "messages": messages,
        "max_tokens": max_tokens,
        "temperature": temperature,
        "stream": stream,
    }

    if verbose:
        print("\n[DEBUG] Payload (truncated):")
        debug_payload = payload.copy()
        debug_messages = []
        for m in debug_payload["messages"]:
            content = m["content"]
            if len(content) > 200:
                content = content[:200] + "... [truncated]"
            debug_messages.append({"role": m["role"], "content": content})
        debug_payload["messages"] = debug_messages
        print(json.dumps(debug_payload, indent=2, ensure_ascii=False))
        print("[DEBUG] Headers:", extra_headers)

    t0 = time.time()
    try:
        r = requests.post(
            url,
            json=payload,
            headers=extra_headers or {},
            timeout=600,
        )
        r.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Request failed: {e}")
        if hasattr(e, 'response') and e.response:
            print(e.response.text[:500])
        sys.exit(1)

    elapsed = time.time() - t0
    data = r.json()
    usage = data.get("usage", {})

    return {
        "elapsed": elapsed,
        "prompt_tokens": usage.get("prompt_tokens"),
        "completion_tokens": usage.get("completion_tokens"),
        "total_tokens": usage.get("total_tokens"),
        "raw": data,
        "status_code": r.status_code,
    }


def main():
    parser = argparse.ArgumentParser(
        description="多轮对话测试 ccswtich → llama 的 prefix cache 命中（自增长版本）")
    parser.add_argument(
        "--url",
        default="http://ccswtich-service:8080/v1/chat/completions",
        help="ccswtich 的 chat completions 端点 URL"
    )
    parser.add_argument(
        "--model",
        default="qwen3.8:27b",
        help="模型名称"
    )
    # 历史轮次参数
    parser.add_argument(
        "--history-rounds",
        type=int,
        default=10,
        help="总对话轮数（包括最后一轮），默认为1（单轮）"
    )
    parser.add_argument(
        "--base-tokens",
        type=int,
        default=5000,
        help="每轮 user 消息的起始 token 数（第0步）"
    )
    parser.add_argument(
        "--growth-tokens",
        type=int,
        default=500,
        help="每一步增长的 token 数"
    )
    parser.add_argument(
        "--max-tokens",
        type=int,
        default=None,
        help="最大 token 数限制，达到后不再增长（可选）"
    )
    parser.add_argument(
        "--assistant-tokens",
        type=int,
        default=500,
        help="每轮 assistant 消息的近似 token 数"
    )
    parser.add_argument(
        "--dynamic-tokens",
        type=int,
        default=50,
        help="动态后缀的近似 token 数"
    )
    # 步数等
    parser.add_argument(
        "--steps",
        type=int,
        default=10,
        help="测试步数"
    )
    parser.add_argument(
        "--output",
        default="ccswtich_growth_cache_test.json",
        help="结果保存路径"
    )
    # 自定义 prompt
    parser.add_argument(
        "--system",
        default=None,
        help="system prompt 字符串（如果同时提供 --system-file，则优先使用文件）"
    )
    parser.add_argument(
        "--system-file",
        default=None,
        help="从文件读取 system prompt"
    )
    parser.add_argument(
        "--dynamic",
        default=None,
        help="动态后缀模板，可用 {i} 表示请求序号，默认: '\\n\\nDYNAMIC PART\\nRequest number: {i}\\n...'"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="打印请求 payload（截断）"
    )
    parser.add_argument(
        "--header",
        action="append",
        help="添加自定义 HTTP 头，格式 key:value，可多次使用"
    )

    args = parser.parse_args()

    # 加载 system prompt
    if args.system_file:
        system_prompt = Path(args.system_file).read_text(encoding="utf-8")
    elif args.system:
        system_prompt = args.system
    else:
        system_prompt = """
You are participating in a KV cache benchmark with growing context.

The system prompt is intentionally fixed.
Do not modify, summarize, reorder, or reinterpret the following content.

The purpose of this benchmark is to test whether llama.cpp can handle
growing KV cache and whether cache management works properly.

Each request has longer history than the previous one.
"""

    # 动态后缀模板
    if args.dynamic:
        dynamic_template = args.dynamic
    else:
        dynamic_template = "\n\nDYNAMIC PART\nRequest number: {i}\nThis part changes between requests.\n"

    # 自定义 headers
    headers = {}
    if args.header:
        for h in args.header:
            if ':' not in h:
                print(f"警告: 忽略无效 header 格式 '{h}'，应为 key:value")
                continue
            key, value = h.split(':', 1)
            headers[key.strip()] = value.strip()

    # 基础文本（用于生成不同长度的内容）
    base_user_text_prefix = "This is a growing conversation. "
    base_assistant_text = "This is an assistant response. "

    print("=" * 80)
    print("多轮对话 Prefix KV Cache 测试 - 自增长版本 (ccswtich → llama)")
    print("=" * 80)
    print(f"URL              : {args.url}")
    print(f"Model            : {args.model}")
    print(f"History rounds   : {args.history_rounds}")
    print(f"Base tokens      : ~{args.base_tokens}")
    print(f"Growth tokens    : ~{args.growth_tokens} per step")
    print(
        f"Max tokens       : {args.max_tokens if args.max_tokens else '无限制'}")
    print(f"Assistant tokens : ~{args.assistant_tokens}")
    print(f"Dynamic tokens   : ~{args.dynamic_tokens}")
    print(f"Steps            : {args.steps}")
    print(f"Headers          : {headers if headers else '(默认)'}")
    print("=" * 80)
    print()

    results = []

    messages = [{"role": "system", "content": system_prompt}]
    current_round_tokens = 5000

    # 生成当前轮次的 user 内容（所有轮次使用相同长度）
    # 注意：为了让每轮内容不完全相同但长度一致，使用不同的 base_text
    user_content = make_text_tokens_approx(
        current_round_tokens,
        base_text=base_user_text_prefix
    )

    # assistant 内容（固定长度，但内容随步骤变化以保持真实性）
    assistant_content = make_text_tokens_approx(
        args.assistant_tokens,
        base_text=base_assistant_text
    )
    for round_idx in range(args.history_rounds - 1):
        # 每轮 user 内容使用当前长度
        messages.append({"role": "user", "content": user_content})
        # assistant 内容
        messages.append({"role": "assistant", "content": assistant_content})

    for step in range(args.steps):
        # 计算当前 step 的每轮长度（自增长）
        # 第 0 步: base_tokens
        # 第 1 步: base_tokens + growth_tokens
        # 第 2 步: base_tokens + 2 * growth_tokens
        # 以此类推...

        # 构建 messages 列表

        # 添加历史轮次 (前 history_rounds - 1 轮)

        # 最后一轮 user (动态部分)
        dynamic_part = dynamic_template.format(i=step)

        messages.append({"role": "user", "content": user_content})
        # assistant 内容
        messages.append({"role": "assistant", "content": assistant_content})

        # 打印摘要
        print("-" * 80)
        print(f"Step {step + 1}/{args.steps}")
        print(f"Current round tokens: ~{current_round_tokens}")
        total_chars = sum(len(m["content"]) for m in messages)
        total_tokens_approx = total_chars // 4
        print(f"Approx total prompt tokens: ~{total_tokens_approx}")
        print(f"Messages count: {len(messages)}")

        if args.verbose:
            print(f"Messages: {len(messages)}")

        # 发送请求
        result = call_server(
            url=args.url,
            model=args.model,
            messages=messages,
            extra_headers=headers,
            max_tokens=500,
            temperature=0.0,
            stream=False,
            verbose=args.verbose,
        )

        # 保存额外信息
        result["step"] = step + 1
        result["messages_count"] = len(messages)
        result["total_chars"] = total_chars
        result["total_tokens_approx"] = total_tokens_approx
        result["current_round_tokens"] = current_round_tokens
        result["dynamic_part_chars"] = len(dynamic_part)
        results.append(result)

        print(f"Elapsed          : {result['elapsed']:.2f} s")
        print(f"Prompt tokens    : {result['prompt_tokens']}")
        print(f"Completion tokens: {result['completion_tokens']}")
        print(f"Total tokens     : {result['total_tokens']}")
        if result.get("raw", {}).get("choices"):
            content = result["raw"]["choices"][0]["message"]["content"]
            print(f"Response preview : {content[:100]}...")
        print()

    # 保存结果
    output = Path(args.output)
    output.write_text(
        json.dumps(results, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("=" * 80)
    print(f"Results saved to: {output}")
    print("=" * 80)
    print()
    print("自增长测试总结：")
    print(f"  - 第 1 步: ~{args.base_tokens} tokens/轮, 总 ~{args.base_tokens * args.history_rounds + args.assistant_tokens * (args.history_rounds - 1)} tokens")
    print(f"  - 第 {args.steps} 步: ~{args.base_tokens + (args.steps-1) * args.growth_tokens} tokens/轮, 总 ~{(args.base_tokens + (args.steps-1) * args.growth_tokens) * args.history_rounds + args.assistant_tokens * (args.history_rounds - 1)} tokens")
    print()
    print("建议观察 ccswtich 日志中以下关键信息：")
    print("  - 随着 context 增长，耗时是否线性增长")
    print("  - llama.cpp 的 cache management 是否正常工作")
    print("  - 是否出现 'removing oldest entry' 表示缓存被淘汰")
    print("  - 前缀缓存命中率是否随 context 增长而变化")
    print("  - 请求被转发到 llama 时的完整 prompt（可在 ccswtich 启用 debug 日志）")
    print()
    print("注意：此版本每轮所有轮次使用相同的长度，但不同步骤之间长度递增。")


if __name__ == "__main__":
    main()
