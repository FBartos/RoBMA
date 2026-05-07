# Workflow

Use this file for repository-level work habits. Tool-specific system rules still
take precedence.

## Planning

- For non-trivial tasks, keep a short visible plan and update it as work changes.
- Do not create repo task files unless the user explicitly asks for durable task tracking.
- Re-plan when new evidence changes the implementation path.
- For complex design choices or unclear maintainer preferences, create a concise `.agents/decisions-<topic>.md` file with context, options, proposed default, and blank `Decision:` fields.

## Scope

- Treat the existing dirty worktree as user work unless proven otherwise.
- Touch only files needed for the request.
- Do not revert unrelated changes.
- Prefer current repo patterns over new abstractions.

## Verification

- Prove changes with the narrowest meaningful command first.
- Use `devtools::test(..., reporter = "llm")` for tests.
- For cached-fit tests, run `devtools::test(filter = "01-", reporter = "llm")` when cache population is needed.
- If tests are not run, state that explicitly and why.

## Instruction Maintenance

- Verify every referenced path/function with `rg` before updating instructions.
- Keep one root `AGENTS.md` plus detailed files under `.agents/instructions/`.
- If maintainer intent is unclear, add a concise item to `.agents/instructions-decisions.md` with a proposed decision.
- After the maintainer adds a `Decision:` line, incorporate it and remove or mark the item resolved.

## Standards

- Find root causes; avoid temporary fixes unless explicitly requested.
- Keep changes small and reversible.
- Surface real ambiguity early, with the concrete tradeoff.
