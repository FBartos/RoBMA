# Instruction Decisions

Use this file for pending maintainer decisions about `AGENTS.md` and files under
`.agents/instructions/`.

Each item should include:

```markdown
## N. Short Decision Title

Context: What is unclear and where it appears.

Why it matters: The tradeoff, risk, or maintenance cost.

Options:
- A: ...
- B: ...

Proposed default: ...

Decision:
```

After `Decision:` is filled in, incorporate the result into the relevant
instruction files, then remove the resolved item or keep a one-line resolved
note.

## Resolved

- 2026-05-05: Treat `engineer-mode.md` as the default engineering prelude for repo work. Because repo files cannot force runtime prompt prepending, `AGENTS.md` now requires reading it at session start and re-reading it after compaction/resume or before non-trivial planning and edits. `engineer-mode.md` now includes the markdown decision-file pattern for complex designs.
- 2026-05-05: Document active p-hacking selection kernels as deferred, hidden from users, and not supported until explicit maintainer approval.
- 2026-05-05: Prefer robust UTF-8 documentation rendering over ASCII-only text. Keep `DESCRIPTION` as `Encoding: UTF-8`, keep vignette `%\VignetteEncoding{UTF-8}`, use UTF-8 names in prose/YAML, and use LaTeX accent escapes in BibTeX.
