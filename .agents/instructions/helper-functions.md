# Helper Function Reuse

Always use existing internal helper functions rather than accessing object internals directly.

## Data Access Helpers

Located in `R/outcome-helpers.R`:

| Helper | Purpose | Use Instead Of |
|--------|---------|----------------|
| `.outcome_data_yi(object)` | Get observed effect-size values on the fitted scale | `object[["data"]][["outcome"]][["yi"]]` |
| `.outcome_data_sei(object)` | Get standard errors | `object[["data"]][["outcome"]][["sei"]]` |
| `.outcome_data_vi(object)` | Get sampling variances | `sei^2` or direct access |
| `.outcome_data_weights(object)` | Get data weights, defaulting to 1 when absent | direct weights access with NULL handling |

## Model Type Helpers

Located in `R/fit.R`:

| Helper | Purpose |
|--------|---------|
| `.is_multilevel(object)` | Check if 3-level model |
| `.is_scale(object)` | Check if location-scale model |
| `.is_weightfunction(object)` | Check if selection model |
| `.is_PET(object)` / `.is_PEESE(object)` | Check PET/PEESE bias models |
| `.is_mods(object)` | Check if location meta-regression is present |
| `.is_weights(object)` | Check if likelihood weights are present |
| `.outcome_type(object)` | Return `"norm"`, `"bin"`, or `"pois"` |
| `.effect_direction(object)` | Return fitted effect direction |

## Evaluation Helpers

Located in `R/evaluate.R`:

| Helper | Purpose |
|--------|---------|
| `.evaluate.brma.tau()` | Extract tau samples (handles multilevel) |
| `.evaluate.brma.mu()` | Extract mu samples |
| `.evaluate.brma.bias_offset()` | Compute PET/PEESE bias offsets |
| `.extract_use_normal()` | Detect posterior rows that can use normal selected-kernel mode |

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
