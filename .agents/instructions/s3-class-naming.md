# S3 Class Naming Convention

For functions that return objects needing print/summary methods, use consistent S3 class naming.

## Pattern

```
{function_name}.{parent_class}
```

## Examples

| Function | Parent Class | S3 Class Name |
|----------|--------------|---------------|
| `summary_heterogeneity()` | `brma` | `summary_heterogeneity.brma` |
| `summary()` | `brma` | `summary.brma` |
| `residuals()` | `brma` | `residuals.brma` |
| `vif()` | `brma` | `vif.brma` |
| `marginal_means()` | `brma` | `marginal_means.brma` |

## Implementation

```r
# Function returns object with class
summary_heterogeneity.brma <- function(object, ...) {

  output <- list(
    estimates = estimates_table
  )

  class(output) <- "summary_heterogeneity.brma"

  return(output)
}

# Print method uses same class
print.summary_heterogeneity.brma <- function(x, ...) {
  cat("\n")
  print(x$estimates)
  cat("\n")
  return(invisible(x))
}
```

## Why This Matters

- Consistent naming makes method dispatch predictable
- Clear relationship between function and its output class
- Follows R convention for S3 generic/method naming
- Makes it easy to find related methods in codebase

## Testing

```r
expect_s3_class(result, "summary_heterogeneity.brma")
```
