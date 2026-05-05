# sampling_bias Parameter Pattern

For RoBMA plotting and analysis functions that use predictions from models with publication bias adjustments, add a `sampling_bias` parameter to control whether bias is incorporated.

## Applies To

- `RoBMA` (ensemble with publication bias)
- `bPET` (PET regression adjustment)
- `bPEESE` (PEESE regression adjustment)
- `bselmodel` (selection/weightfunction models)

## Parameter Definition

```r
#' @param sampling_bias whether publication bias should be incorporated into the
#' predicted effect. Defaults to \code{TRUE}. When \code{TRUE} and the model
#' includes PET/PEESE, incorporates the expected bias from these regression
#' adjustments into predictions. When \code{FALSE}, shows bias-adjusted
#' (corrected) predictions.
```

## Implementation

The `sampling_bias` parameter maps to `predict.brma()`'s `bias_adjusted` with **inverted logic**:

```r
# In function signature
sampling_bias = TRUE

# When calling predict.brma()
pred_samples <- predict.brma(
  object        = x,
  newdata       = newdata,
  type          = "terms",
  bias_adjusted = !sampling_bias,  # inverted!
  quiet         = TRUE
)
```

## Behavior

| `sampling_bias` | `bias_adjusted` | Effect |
|-----------------|-----------------|--------|
| `TRUE` (default) | `FALSE` | PET/PEESE terms ARE added → shows biased predictions (what you'd observe) |
| `FALSE` | `TRUE` | PET/PEESE terms NOT added → shows corrected predictions |

## Rationale

- Default `TRUE` matches user expectation: "include the bias in what I see"
- Consistent with `funnel.brma()` which has `sampling_bias` parameter
- The inversion handles the semantic difference between "include bias" (user-facing) vs "adjust for bias" (internal)

## Example Usage

```r
# Show regression line with publication bias effect included (default)
regplot(fit, sampling_bias = TRUE)

# Show bias-corrected regression line
regplot(fit, sampling_bias = FALSE)
```
