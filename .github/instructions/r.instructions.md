---
description: 'R language and document formats (R, Rmd, Quarto): coding standards and Copilot guidance for idiomatic, safe, and consistent code generation.'
applyTo: '**/*.R, **/*.r, **/*.Rmd, **/*.rmd, **/*.qmd'
---

# R Programming Language Instructions

## Purpose

Help GitHub Copilot generate idiomatic, safe, and maintainable R code across projects.

## Core Conventions

- **Match the project’s style.** Follow the style in the project.
- **Prefer clear, vectorized code.** Keep functions small and avoid hidden side effects.
- **Qualify non-base functions in examples/snippets**, e.g., `dplyr::mutate()`, `stringr::str_detect()`.
- **Naming:** `lower_snake_case` for objects/files; use dots to dispatch different function types (and in S3 classes).
- **Side effects:** Never call `setwd()`; prefer project-relative paths (e.g., `here::here()`).
- **Validation:** Validate and constrain user inputs; use the predefined `check_bool()`, `check_char()`, `check_real()` ... functions.

### Pipe Operators

- **Never use pipe:** Always assign values using an arror `<-`

## Performance Considerations

- **Profiling:** Use `profvis::profvis()` to identify performance bottlenecks in your code. Profile before optimizing.
- **Caching:** Use `memoise::memoise()` to cache expensive function results. Particularly useful for repeated API calls or complex computations.
- **Vectorization:** Prefer vectorized operations over loops. Use `apply()` family for remaining iteration needs.

## Tooling & Quality

- **Pre-commit:** consider `precommit` hooks to lint/format automatically.
- **Docs:** roxygen2 for exported functions (`@param`, `@return`, `@examples`).
- **Tests:** prefer small, pure, composable functions that are easy to unit test.

## Data Wrangling & I/O

- **Data frames:** Use base `data.frame()`
- **Iteration:** Prefer type-stable, vectorized patterns such as `vapply()` (for atomic outputs). Use `for` loops when when they improve clarity or performance.
- **Strings & Dates:** Use clear base helpers (e.g., `nchar()`, `substr()`, `as.Date()` with explicit format).
- **I/O:** prefer explicit, typed readers (e.g., `readr::read_csv()`); make parsing assumptions explicit.

## Error Handling

- Use `stop(..., .call = FALSE)` / `warning()`.
- For recoverable operations:
- Use `tryCatch()` in base R for fine-grained control.

## Security Best Practices

- **Command execution:** Prefer `processx::run()` or `sys::exec_wait()` over `system()`; validate and sanitize all arguments.
- **File paths:** Normalize and sanitize user-provided paths (e.g., `fs::path_sanitize()`), and validate against allowlists.
- **Credentials:** Never hardcode secrets. Use env vars (`Sys.getenv()`), config outside VCS, or `keyring`.

## Copilot-Specific Guidance

- Suggest vectorized solutions over loops when idiomatic.
- Prefer small helper functions over long pipelines.
- When multiple approaches are equivalent, prefer readability and type stability and explain the trade-offs.

---

## Minimal Examples

```r
scores <- data.frame(id = 1:5, x = c(1, 3, 2, 5, 4))
safe_log <- function(x) tryCatch(log(x), error = function(e) NA_real_)
scores$z <- vapply(scores$x, safe_log, numeric(1))

# Example reusable helper with roxygen2 doc
#' Compute the z-score of a numeric vector
#' @param x A numeric vector
#' @return Numeric vector of z-scores
#' @examples z_score(c(1, 2, 3))
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
```
