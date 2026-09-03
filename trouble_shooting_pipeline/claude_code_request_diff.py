#!/usr/bin/env python3
"""
Locate where a llama.cpp KV/prompt cache is invalidated between two Claude Code
request bodies.

llama.cpp's cache is a *token* prefix cache: the KV cache is built token by
token in order, and the reusable prefix is the longest run of matching tokens
(the modern `llama_kv_cache` uses a hash chain over the token sequence). So
"where does the cache invalidate" == the first index where the two *tokenized
prompts* diverge.

Method
------
1. Render each request into the prompt text the backend would tokenize
   (system -> tools -> messages, in wire order). Rendering is deterministic
   and identical on both sides, so the divergence boundary is valid.
2. Tokenize with the model tokenizer (--tokenizer) and take the token-level
   longest common prefix. This IS the cached prefix length (in tokens).
3. Map the divergent token back to a JSON path (messages[i].content[j]) via a
   char-offset -> span table captured during rendering.

If no tokenizer is given, a char-level LCP on the rendered text is used as a
proxy (chars track tokens closely). The *boundary location* is robust either
way; the *exact cached-token count* is only exact when the tokenizer matches
the model.
"""

import argparse
import bisect
import json
import sys
from typing import Any


def load_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Request extraction
# ---------------------------------------------------------------------------

def get_request_body(obj: Any) -> Any:
    """
    Return the request body (the object carrying system/tools/messages) from
    several common dump formats.

    1. Pure request JSON:   {"model": ..., "messages": ...}
    2. Wrapped / mitmproxy: {"request": {"content": ...}} where `content` is
       a raw JSON string or a dict.
    """

    if isinstance(obj, dict):
        if any(key in obj for key in
               ["messages", "system", "tools", "model", "max_tokens"]):
            return obj
        if "request" in obj:
            req = obj["request"]
            if isinstance(req, dict):
                content = req.get("content")
                if content is not None:
                    if isinstance(content, str):
                        try:
                            return json.loads(content)
                        except json.JSONDecodeError:
                            return content
                    return content
                return req
    return obj


def as_list(content: Any) -> list:
    """Normalize a message `content` field to a list of blocks."""
    if isinstance(content, str):
        return [{"type": "text", "text": content}]
    if isinstance(content, list):
        return content
    return [content]


def content_block_text(block: Any) -> str:
    """
    A deterministic text form of one content block — the cache-relevant
    content. text/thinking blocks yield their text; tool_result yields its
    raw content; other blocks (tool_use) are serialized as-is so any change
    to a tool call is still detected.
    """

    if isinstance(block, dict):
        btype = block.get("type", "")

        if btype == "text":
            return block.get("text", "")
        if btype in ("thinking", "redacted_thinking"):
            return block.get("thinking", "") + block.get("data", "")
        if btype == "tool_result":
            c = block.get("content")
            if isinstance(c, str):
                return c
            if isinstance(c, list):
                return "".join(content_block_text(x) for x in c)
            return json.dumps(c, ensure_ascii=False, sort_keys=False,
                              separators=(",", ":"))
        return json.dumps(block, ensure_ascii=False, sort_keys=False,
                          separators=(",", ":"))

    return str(block)


def canonical_block(block: Any) -> str:
    """Order-preserving serialization for a whole field (e.g. tools)."""
    return json.dumps(block, ensure_ascii=False, sort_keys=False,
                      separators=(",", ":"))


# ---------------------------------------------------------------------------
# Rendering: request -> prompt text + span table
# ---------------------------------------------------------------------------

def render_prompt(req: dict):
    """
    Render the request into a single prompt string in wire order, and return
    (text, spans) where spans is a list of (start, end, json_path) covering
    the text. Any char offset can be mapped back to a JSON path via the
    spans.

    Order is system -> tools -> messages. All of this is deterministic and
    identical on both requests, so it does not move the divergence boundary.
    """

    parts = []  # (text, path)

    def seg(text: str, path: str):
        parts.append((text, path))

    system = req.get("system")
    if isinstance(system, str):
        seg(system, "$.system")
    elif isinstance(system, list):
        for j, b in enumerate(system):
            seg(content_block_text(b), f"$.system[{j}]")

    tools = req.get("tools")
    if tools:
        seg(canonical_block(tools), "$.tools")

    for i, m in enumerate(req.get("messages", [])):
        role = m.get("role", "?") if isinstance(m, dict) else "?"
        seg(f"\n<<MSG {i} role={role}>>\n", f"$.messages[{i}].role")
        for j, blk in enumerate(as_list(m.get("content"))):
            seg(content_block_text(blk), f"$.messages[{i}].content[{j}]")

    text = "".join(t for t, _ in parts)
    spans = []
    pos = 0
    for t, path in parts:
        spans.append((pos, pos + len(t), path))
        pos += len(t)

    return text, spans


def locate_span(offset: int, spans: list, starts: list) -> str:
    """Map a char offset to the JSON path of the span containing it."""
    idx = bisect.bisect_right(starts, offset) - 1
    if idx < 0:
        return spans[0][2] if spans else "?"
    return spans[idx][2]


# ---------------------------------------------------------------------------
# LCP helpers
# ---------------------------------------------------------------------------

def lcp_chars(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def hr(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def quick_prefix_status(req_a: dict, req_b: dict) -> None:
    """Cheap cross-ref: is the system / tools base identical?"""
    hr("0. BASE CROSS-REF  (system / tools)")
    for field in ("system", "tools"):
        fa = req_a.get(field)
        fb = req_b.get(field)
        if canonical_block(fa) == canonical_block(fb):
            print(f"  {field}: identical")
        else:
            print(f"  {field}: *** DIFFERS ***  (this is an early cache-break)")


# ---------------------------------------------------------------------------
# Primary: token-level boundary
# ---------------------------------------------------------------------------

def token_boundary(text_a, text_b, spans_a, spans_b, tokenizer) -> bool:
    """
    Token-level LCP over the rendered prompts. Returns True if a token
    boundary was fully resolved (so callers can skip the char fallback).

    On success prints the cached token count, the divergent token position,
    decoded context, and the mapped JSON path.
    """

    try:
        ea = tokenizer(text_a, return_offsets_mapping=True,
                       add_special_tokens=False)
        eb = tokenizer(text_b, return_offsets_mapping=True,
                       add_special_tokens=False)
        ta, tb = ea["input_ids"], eb["input_ids"]
        oa, ob = ea["offset_mapping"], eb["offset_mapping"]
    except Exception as e:
        print(f"  [token-level unavailable: {e} — falling back to char-level]")
        return False

    n = min(len(ta), len(tb))
    i = 0
    while i < n and ta[i] == tb[i]:
        i += 1

    starts_a = [s for s, _, _ in spans_a]
    starts_b = [s for s, _, _ in spans_b]

    print(f"A tokens : {len(ta):,}")
    print(f"B tokens : {len(tb):,}")
    print(f"Token LCP: {i:,}  (cached prefix ~ {100 * i / max(len(tb), 1):.1f}% of B)")

    if i == n and len(ta) == len(tb):
        print("\n> Prompts identical -> full cache hit.")
        return True

    if i == n:
        longer = "B" if len(tb) > len(ta) else "A"
        print(f"\n> One prompt is a strict prefix of the other ({longer} longer)"
              f" -> pure append, full cache hit for the shared prefix.")
        return True

    # First divergent token is index i; tokens 0..i-1 are the cached prefix.
    off_b = ob[i][0]
    off_a = oa[i][0]
    path_b = locate_span(off_b, spans_b, starts_b)
    path_a = locate_span(off_a, spans_a, starts_a)

    ctx = 20
    s0 = max(0, i - ctx)
    e_a = min(len(ta), i + ctx)
    e_b = min(len(tb), i + ctx)

    print(f"\nFirst divergent token = {i:,}")
    print(f"  A char offset = {off_a:,}  ->  {path_a}")
    print(f"  B char offset = {off_b:,}  ->  {path_b}")

    print(f"\n--- A context (tokens {s0}..{e_a - 1}) ---")
    print(repr(tokenizer.decode(ta[s0:e_a], skip_special_tokens=False)))
    print(f"\n--- B context (tokens {s0}..{e_b - 1}) ---")
    print(repr(tokenizer.decode(tb[s0:e_b], skip_special_tokens=False)))

    print(f"\n> CONCLUSION: cache invalidated at token {i:,} -> {path_b}")
    return True


def char_boundary(text_a, text_b, spans_b) -> None:
    """Fallback: char-level LCP on the rendered text, mapped to a JSON path."""

    pos = lcp_chars(text_a, text_b)
    starts_b = [s for s, _, _ in spans_b]
    path = locate_span(pos, spans_b, starts_b)

    print(f"A chars : {len(text_a):,}")
    print(f"B chars : {len(text_b):,}")
    print(f"Char LCP: {pos:,}  (~ {100 * pos / max(len(text_b), 1):.1f}% of B)")

    if pos == min(len(text_a), len(text_b)) and len(text_a) == len(text_b):
        print("\n> Prompts identical -> full cache hit.")
        return
    if pos == min(len(text_a), len(text_b)):
        longer = "B" if len(text_b) > len(text_a) else "A"
        print(f"\n> One prompt is a strict prefix of the other ({longer} longer)"
              f" -> pure append, full cache hit for the shared prefix.")
        return

    ctx = 300
    s0 = max(0, pos - ctx)
    print(f"\nFirst divergence at char {pos:,}  ->  {path}")
    print(f"\n--- A around pos {pos:,} ---")
    print(repr(text_a[s0:pos + ctx]))
    print(f"\n--- B around pos {pos:,} ---")
    print(repr(text_b[s0:pos + ctx]))

    if pos < min(len(text_a), len(text_b)):
        print(f"\n  A[{pos}] = {text_a[pos]!r}")
        print(f"  B[{pos}] = {text_b[pos]!r}")
    print(f"\n> CONCLUSION (char proxy): cache invalidated at char {pos:,} -> {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Locate the llama.cpp cache-invalidation boundary between "
                    "two Claude Code request bodies (token-level, with JSON "
                    "path mapping).")

    parser.add_argument("request_a", help="First request JSON")
    parser.add_argument("request_b", help="Second request JSON")

    parser.add_argument("--tokenizer", default=None,
                        help="HuggingFace tokenizer matching the llama.cpp "
                             "model, e.g. Qwen/Qwen3-30B-A3B. Without it, a "
                             "char-level proxy is used.")

    args = parser.parse_args()

    print("=" * 80)
    print("Loading requests")
    print("=" * 80)

    obj_a = load_json(args.request_a)
    obj_b = load_json(args.request_b)

    req_a = get_request_body(obj_a)
    req_b = get_request_body(obj_b)

    print(f"A: {args.request_a}")
    print(f"B: {args.request_b}")

    if not (isinstance(req_a, dict) and isinstance(req_b, dict)):
        print("\n[ERROR] Neither input looks like a request body "
              "(no messages/system/tools). Aborting.")
        sys.exit(1)

    quick_prefix_status(req_a, req_b)

    text_a, spans_a = render_prompt(req_a)
    text_b, spans_b = render_prompt(req_b)

    if args.tokenizer:
        try:
            from transformers import AutoTokenizer
        except ImportError:
            print("\n[WARN] transformers not installed. "
                  "Install with:  pip install transformers")
            AutoTokenizer = None

        if AutoTokenizer is not None:
            hr("1. TOKEN-LEVEL CACHE BOUNDARY")
            tokenizer = AutoTokenizer.from_pretrained(
                args.tokenizer, trust_remote_code=True)
            resolved = token_boundary(text_a, text_b, spans_a, spans_b,
                                      tokenizer)
            if not resolved:
                hr("2. CHAR-LEVEL FALLBACK")
                char_boundary(text_a, text_b, spans_b)
        else:
            hr("2. CHAR-LEVEL (no tokenizer)")
            char_boundary(text_a, text_b, spans_b)
    else:
        print("\n[note] --tokenizer not given; using char-level proxy. "
              "For an exact token count, pass --tokenizer <model>.")
        hr("CHAR-LEVEL CACHE BOUNDARY (proxy)")
        char_boundary(text_a, text_b, spans_b)


if __name__ == "__main__":
    main()
