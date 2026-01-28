# R Development Instructions

Reference this file when working with R code (`.R`, `.Rmd`, `.qmd` files).

## Core Conventions

- **Match the project's style.** Follow the style in the project.
- **Prefer clear, vectorized code.** Keep functions small and avoid hidden side effects.
- **Qualify non-base functions**, e.g., `BayesTools::check_bool()`, `stats::sd()`.
- **Naming:** `lower_snake_case` for objects/files; use dots for S3 method dispatch.
- **Side effects:** Never call `setwd()`; prefer project-relative paths.
- **Validation:** Use the predefined `BayesTools::check_bool()`, `check_char()`, `check_real()` functions.

### Assignment & Pipes

- **Always use `<-`** for assignment, never `=`
- **Never use pipes** (`|>` or `%>%`); always assign intermediate values

### Booleans

- **Always use `TRUE`/`FALSE`**, never `T`/`F`

## Performance Considerations

- **Profiling:** Use `profvis::profvis()` to identify performance bottlenecks.
- **Caching:** Use `memoise::memoise()` for expensive function results.
- **Vectorization:** Prefer vectorized operations over loops. Use `vapply()` for type-stable iteration.

## Tooling & Quality

- **Docs:** roxygen2 for exported functions (`@param`, `@return`, `@examples`).
- **Tests:** Prefer small, pure, composable functions that are easy to unit test.

## Data Wrangling & I/O

- **Data frames:** Use base `data.frame()`
- **Iteration:** Prefer type-stable, vectorized patterns such as `vapply()` (for atomic outputs). Use `for` loops when they improve clarity or performance.
- **Strings & Dates:** Use clear base helpers (e.g., `nchar()`, `substr()`, `as.Date()` with explicit format).

## Error Handling

- Use `stop(..., call. = FALSE)` / `warning()`.
- Use `tryCatch()` for recoverable operations.

## Security Best Practices

- **Command execution:** Prefer `processx::run()` over `system()`; validate all arguments.
- **File paths:** Normalize and sanitize user-provided paths.
- **Credentials:** Never hardcode secrets. Use env vars (`Sys.getenv()`).
