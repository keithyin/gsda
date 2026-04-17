# Terminal-Style Output Display for CLI Tools

## Problem
When viewing CLI tool outputs (especially Polars DataFrame tables with ASCII formatting) in the ResultViewer, the displayed results didn't match terminal output. Table column alignment was broken.

## Root Causes

1. **CSS `white-space: pre-wrap`** caused browser to wrap long lines, breaking table alignment
2. **`v-html` directive with `highlightSearch()`** processed text as HTML, interfering with raw text display
3. **Extra line numbers** added to each line disturbed the original output formatting
4. **`UTF8` table formatting** instead of `ASCII` for Polars output

## Solution

### 1. Terminal-Style CSS
```css
.terminal-output {
  font-family: 'Courier New', 'Courier', monospace;
  background: #1e1e1e;
  color: #d4d4d4;
  white-space: pre;        /* Preserve all whitespace */
  word-wrap: normal;       /* Don't wrap lines */
  overflow-x: auto;        /* Horizontal scroll */
  overflow-y: auto;
  display: block;
}
```

### 2. Use Pure Text Rendering
```vue
<pre class="output-text terminal-output">
  {{ result.stdout }}
</pre>
```
- Remove `v-html` directive
- Remove line number display
- Use `{{ }}` for text interpolation

### 3. CLI Tool ASCII Formatting
```python
def pretty_print(obj):
    pl.Config.set_tbl_formatting("ASCII")  # Not "UTF8"
    width = shutil.get_terminal_size().columns
    pl.Config.set_tbl_width(width)
    print(obj)
```

## Files Modified
- `src/gseda/server/frontend/src/components/ResultViewer.vue` - Output display component
- `src/gseda/bam_ana/bam_basic_stat.py` - Table formatting configuration

## Key Takeaways
- Always use `white-space: pre` for preserving code/table formatting in browser
- Avoid `v-html` for displaying raw terminal output
- Strip extra metadata (like line numbers) from raw output display
- Use ASCII formatting for terminal-compatible output
