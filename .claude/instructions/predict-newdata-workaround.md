# predict.brma Newdata Workaround

When calling `predict.brma()` with `newdata` for `type="terms"` or `type="estimate"` predictions (which only need moderator values for formula evaluation), dummy outcome columns must still be provided because `.prepare_newdata()` validates their presence.

## Why This Happens

The `.prepare_newdata()` function requires outcome columns (yi/sei for normal, ai/ci/n1i/n2i for binomial, etc.) even though they're not used for `type="terms"` predictions. This is because the same validation code path is used for all prediction types.

## Solution

Add dummy outcome values with zeros. The `skip_validation = TRUE` flag in `.prepare_newdata()` allows zero values:

```r
# Determine outcome type
outcome_type <- .outcome_type(x)

if (outcome_type == "norm") {
  # Normal models: add dummy yi/sei
  newdata[["yi"]]  <- rep(0, n_pred)
  newdata[["sei"]] <- rep(0, n_pred)
} else if (outcome_type == "bin") {
  # Binomial GLMM: add dummy counts (zeros OK with skip_validation)
  newdata[["ai"]]  <- rep(0, n_pred)
  newdata[["ci"]]  <- rep(0, n_pred)
  newdata[["n1i"]] <- rep(0, n_pred)
  newdata[["n2i"]] <- rep(0, n_pred)
} else if (outcome_type == "pois") {
  # Poisson GLMM: add dummy counts
  newdata[["x1i"]] <- rep(0, n_pred)
  newdata[["x2i"]] <- rep(0, n_pred)
  newdata[["t1i"]] <- rep(0, n_pred)
  newdata[["t2i"]] <- rep(0, n_pred)
}
```

## Important Notes

1. These dummy values are **not used** for `type="terms"` or `type="estimate"` predictions
2. The `skip_validation = TRUE` flag in `.prepare_newdata()` allows zero/invalid values
3. For `type="response"` predictions, real outcome data is required

## When to Use

- `regplot()` - needs predictions at new moderator values
- Any function generating predictions for visualization at arbitrary moderator values
- NOT needed when using `newdata = NULL` (existing data) or `newdata = TRUE` (aggregated)
