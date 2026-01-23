---
description: Persistent instructions and preferences for the AI agent
---

# Agent Instructions

## Code Style Preferences

### Assignment Alignment
- **CRITICAL**: Always follow the assignment alignment rules defined in `/r-code-style`.
- When modifying blocks of code, ensure that `<-` and `=` operators are aligned with existing code in that block.
- Example of preferred assignment alignment:
  ```r
  fit_metafor <- info[[name]][["metafor"]]
  fit_brma    <- fits[[name]]
  ```
- Do NOT remove existing alignment even if it seems redundant.

### Assignment Operators
- Use `<-` for assignment, NOT `=`.

### Spacing
- Maintain 2-space indentation for all R code.
- Always include an empty line at the start of a function body (immediately after the opening brace `{`).

## Debugging Workflow
- When tests fail, check for namespace issues, especially for generics that might be masked by or masking other packages (e.g., `loo`).
