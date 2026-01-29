# Metafor Reference Comparison Testing

When writing tests that compare RoBMA/brma results against metafor:

## Pattern

1. **Save metafor fit during model creation** (in `test-01-*.R` files):
   ```r
   fit_metafor <- metafor::rma(yi, vi, data = dat)
   save_fit("model_name", fit_brma, info = list(metafor = fit_metafor))
   ```

2. **Load both fits in comparison tests** (in `test-02-*.R` files):
   ```r
   fits <- lapply(list_fits(), load_fit)
   info <- lapply(list_fits(), load_info)
   names(fits) <- list_fits()
   names(info) <- list_fits()

   fit_brma    <- fits[[name]]
   fit_metafor <- info[[name]][["metafor"]]
   ```

3. **Compare explicitly against metafor output**:
   ```r
   expect_equal(brma_result$estimates["tau", "Mean"],
                sqrt(fit_metafor$tau2), tolerance = 0.05,
                info = "brma tau should match metafor")
   ```

## Why This Matters

- Tests validate statistical correctness against established reference implementation
- Cached fits avoid recomputation during test iteration
- Info objects store any additional comparison data needed

## Files Using This Pattern

- `tests/testthat/test-01-*.R` - Generate and save fits with metafor info
- `tests/testthat/test-02-*.R` - Load fits and compare against metafor
- `tests/testthat/common-functions.R` - `save_fit()`, `load_fit()`, `load_info()` helpers
