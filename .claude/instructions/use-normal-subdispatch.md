# use_normal Subdispatch Pattern for Performance Optimization

When weighted normal functions need to handle mixed samples (some with omega = all 1s, 
some with actual weights), use the `use_normal` subdispatch pattern for efficiency.

## Pattern

The `use_normal` parameter is a logical vector of length S (number of posterior samples)
indicating which samples can use the fast normal path:

- `use_normal[i] = TRUE` → use `dnorm()`/`pnorm()`/`rnorm()` (fast path)
- `use_normal[i] = FALSE` → use full weighted computation (slow path)

## Implementation Template

```r
.dwnorm_fast.ss.matrix <- function(x, mean, sd, omega, crit_x, log = FALSE,
                                    use_normal = NULL) {
  S <- length(mean)
  
  # FAST PATH: If use_normal provided, subdispatch per-row
  if (!is.null(use_normal)) {
    
    log_lik <- rep(NA_real_, S)
    
    # Normal path (fast)
    if (any(use_normal)) {
      log_lik[use_normal] <- stats::dnorm(x, mean = mean[use_normal],
                                          sd = sd[use_normal], log = TRUE)
    }
    
    # Weighted path (slow - recursive call with use_normal = NULL)
    if (any(!use_normal)) {
      idx <- which(!use_normal)
      log_lik[idx] <- .dwnorm_fast.ss.matrix(
        x = x, mean = mean[idx], sd = sd[idx],
        omega = omega[idx, , drop = FALSE], crit_x = crit_x,
        log = TRUE, use_normal = NULL  # <-- Triggers full computation
      )
    }
    
    if (log) return(log_lik) else return(exp(log_lik))
  }
  
  # NO use_normal provided: full weighted computation for ALL rows
  # ... existing code unchanged ...
}
```

## Key Points

- **Default `NULL`**: If `use_normal` not provided, all rows use weighted computation
- **Recursive call**: Avoids code duplication by calling self with `use_normal = NULL`
- **Matrix subsetting**: Use `omega[idx, , drop = FALSE]` to maintain matrix structure
- **Correct positioning**: Results placed back at correct indices

## Where Applied

| Function | Fast path | Slow path |
|----------|-----------|-----------|
| `.dwnorm_fast.ss.matrix()` | `dnorm()` | full weighted density |
| `.pwnorm_fast.ss()` | `pnorm()` | full weighted CDF |
| `.rwnorm_fast.ss()` | `rnorm()` | rejection sampling |
| `.rwnorm_true_fast.ss()` | `rnorm(mean, tau)` | rejection sampling |

## Computing `use_normal`

Use `.extract_use_normal()` in `evaluate.R` to compute this from bias_indicator:

```r
use_normal <- .extract_use_normal(object)
```

## Files

- `R/distributions.R` - Core functions with `use_normal` parameter
- `R/evaluate.R` - `.extract_use_normal()` helper
- `R/pdf.R`, `R/cdf.R`, `R/rng.R` - Wrapper functions passing `use_normal` through
