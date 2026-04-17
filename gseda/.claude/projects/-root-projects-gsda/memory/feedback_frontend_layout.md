# Frontend Flexbox Layout Best Practices for Vue Components

**Rule:** When creating Vue component containers with `display: flex; flex-direction: column;`, always use `min-height: 0` and `flex: 1` instead of `height: 100%` to prevent overflow issues.

**Why:** The combination of `height: 100%` with flexbox column layout can cause:
1. Child elements (like form buttons) to be hidden or cut off
2. Container to expand beyond viewport
3. z-index stacking issues where result viewers or headers overlap form content
4. Scrollbar not appearing for long content

**How to apply:**

## ✅ Recommended Pattern (Correct)
```css
.form-container {
  overflow-y: auto;        /* Enable scrolling for long content */
  display: flex;
  flex-direction: column;
  min-height: 0;           /* Critical: allows flex items to shrink properly */
  flex: 1;                 /* Allow container to grow/fill available space */
  max-height: calc(100vh - 400px);  /* Reserve space for header, results */
  z-index: 1;              /* Layer control when overlays exist */
}
```

## ❌ Problematic Pattern (Avoid)
```css
.form-container {
  height: 100%;            /* Causes flex shrinking issues */
  overflow-y: auto;
  display: flex;
  flex-direction: column;
}
```

## Key Principles

1. **min-height: 0** - Essential for flex items containing other flex items or content with intrinsic size
2. **flex: 1** - Allows container to fill available vertical space while respecting max-height
3. **max-height: calc(100vh - X)** - Reserve space for other UI elements (headers, result viewers, etc.)
4. **z-index** - Use when components may overlap (result viewers, dialogs, etc.)

## Form Button Visibility Pattern

For form action buttons (like "Submit" or "Execute"), ensure they stay visible:

```css
.form-actions {
  margin-top: 30px;
  padding-top: 20px;
  border-top: 1px solid #e4e7ed;
  flex-shrink: 0;  /* Prevent button from being squeezed */
}
```

```typescript
// In component script, add scroll helper:
const scrollToActions = () => {
  const formContainer = document.querySelector('.form-container') as HTMLElement
  if (formContainer) {
    const actions = formContainer.querySelector('.form-actions') as HTMLElement
    if (actions) {
      actions.scrollIntoView({ behavior: 'smooth', block: 'center' })
    }
  }
}

// Call before form submission:
const onSubmit = async () => {
  await nextTick()
  scrollToActions()
  // ... rest of submission logic
}
```

## Component Hierarchy Best Practices

When a parent component contains multiple sections (e.g., form + result viewer):

```css
/* Parent panel */
.tool-panel {
  flex: 1;
  display: flex;
  flex-direction: column;
  overflow: hidden;  /* Prevent scroll at parent level */
}

/* Form section */
.form-container {
  flex: 1;
  min-height: 0;
  max-height: calc(100vh - 400px);
  z-index: 1;
}

/* Result section */
.result-viewer {
  flex-shrink: 0;   /* Prevent result from shrinking */
  z-index: 10;      /* Higher than form when overlaying */
  position: relative;
  margin-top: 20px;
  margin-bottom: 20px;
}
```

## Common Symptoms of Flex Layout Issues

| Symptom | Root Cause | Solution |
|---------|------------|----------|
| Form buttons not visible | Container overflows viewport | Add `max-height` |
| Submit button cut off | No scroll on container | Add `overflow-y: auto` |
| Buttons get squished | `flex: none` or missing `flex-shrink: 0` | Add `flex-shrink: 0` to button container |
| Result viewer overlaps form | Missing z-index layering | Add `z-index` to both components |
| Form content not scrollable | Missing flex properties | Add `min-height: 0; flex: 1` |

## Testing Checklist

- [ ] All form elements visible without scrolling past viewport
- [ ] Submit/button always accessible (scroll to it if needed)
- [ ] Result viewer doesn't hide form inputs
- [ ] Scrollbar appears when content exceeds max-height
- [ ] No overflow causing horizontal scrollbars
- [ ] z-index prevents unintended overlaps
