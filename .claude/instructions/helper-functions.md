# Helper Function Reuse

Always use existing internal helper functions rather than accessing object internals directly.

## Data Access Helpers

Located in `R/brma.residuals.R`:

| Helper | Purpose | Use Instead Of |
|--------|---------|----------------|
| `.outcome_data_sei(object)` | Get standard errors | `object[["data"]][["outcome"]][["sei"]]` |
| `.outcome_data_vi(object)` | Get sampling variances | `sei^2` or direct access |

## Model Type Helpers

Located in various `R/brma.*.R` files:

| Helper | Purpose |
|--------|---------|
| `.is_multilevel(object)` | Check if 3-level model |
| `.is_scale(object)` | Check if location-scale model |
| `.is_weightfunction(object)` | Check if selection model |

## Evaluation Helpers

Located in `R/brma.evaluate.R`:

| Helper | Purpose |
|--------|---------|
| `.evaluate.brma.tau()` | Extract tau samples (handles multilevel) |
| `.evaluate.brma.mu()` | Extract mu samples |

## Why This Matters

- Helpers handle edge cases (different outcome types, missing data)
- Consistent behavior across all model types
- Easier maintenance when data structures change
- Avoids duplicating logic across functions

## Example

```r
# Good: Use helper
vi <- .outcome_data_vi(object)

# Bad: Direct access
vi <- object[["data"]][["outcome"]][["sei"]]^2
```
