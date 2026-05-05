# BayesTools forward API update for RoBMA

Date: 2026-05-04

This is the implementation follow-up to `docs/bayestools-forward-spec.md`.
The requested BayesTools APIs have been added locally in `C:\R-Packages\BayesTools`.
RoBMA can now stop calling the relevant `BayesTools:::` private helpers.

## BayesTools APIs now available

### Formula intercept repair

Use:

```r
BayesTools::formula_add_intercept(formula)
```

Behavior:

- Exported public API.
- Accepts and returns a formula.
- Preserves the formula environment.
- Removes top-level no-intercept encodings such as `- 1`, `+ 0`, `0 +`, and `x + -1`.
- Does not alter transformed or special calls such as `I(x - 1)` or `offset(z - 1)`.

RoBMA connection:

- In `R/fit.R::.create_fit_formula_list()`, replace:

```r
formula <- BayesTools:::.add_intercept_to_formula(formula)
```

with:

```r
formula <- BayesTools::formula_add_intercept(formula)
```

Keep RoBMA's warning and keep setting:

```r
attr(formula, "log(intercept)") <- TRUE
```

after the intercept repair for scale models.

Suggested tests:

- Scale formula `~ x - 1`.
- Scale formula `~ 0 + x`.
- Formula environment preservation.
- `I(x - 1)` and `offset(z - 1)` remain unchanged.

### Transformed prior plotting

Use the general wrapper:

```r
BayesTools::plot_transformed_prior(
  prior_list,
  column_names,
  formula_scale = NULL,
  parameter,
  n_points = 1000,
  x_range = NULL,
  transformation = NULL,
  transformation_arguments = NULL,
  transformation_settings = FALSE,
  plot_type = c("base", "ggplot"),
  par_name = NULL,
  ...
)
```

Behavior:

- Uses BayesTools' internal transformed-prior density machinery.
- Supports continuous densities and point-mass/spike priors.
- Returns base plot metadata invisibly for `plot_type = "base"` and a `ggplot` object for `plot_type = "ggplot"`.
- Returns `NULL` when the requested parameter is an identity transform and no output transformation is requested.
- Accepts `formula_scale` in the same nested shape as `attr(fit, "formula_scale")`.

RoBMA connection:

- In `R/plot.R::.plot_prior_unstandardized()`, replace the three private calls:

```r
BayesTools:::.generate_transformed_prior_densities(...)
BayesTools:::.prior_linear_density_to_plot_data(...)
BayesTools:::.plot_prior_list.both(...)
```

with one call to `BayesTools::plot_transformed_prior()`.

Suggested implementation shape:

```r
par_name <- dots[["par_name"]]
n_points <- if (!is.null(dots[["n_points"]])) dots[["n_points"]] else 1000
BayesTools::check_int(n_points, "n_points", lower = 2)

plot_dots <- dots
plot_dots[c(
  "par_name", "n_points", "n_samples", "force_samples", "x_range_quant",
  "individual", "show_figures", "rescale_x"
)] <- NULL

plot <- do.call(
  BayesTools::plot_transformed_prior,
  c(
    list(
      prior_list              = formula_info[["prior_list"]],
      column_names            = formula_info[["column_names"]],
      formula_scale           = formula_info[["formula_scale"]],
      parameter               = parameter,
      n_points                = n_points,
      x_range                 = dots[["xlim"]],
      transformation          = dots[["transformation"]],
      transformation_arguments = dots[["transformation_arguments"]],
      transformation_settings = if (!is.null(dots[["transformation_settings"]])) {
        dots[["transformation_settings"]]
      } else {
        FALSE
      },
      plot_type               = plot_type,
      par_name                = par_name
    ),
    plot_dots
  )
)

if (is.null(plot)) {
  return(NULL)
}
return(plot)
```

After this change, `.plot_prior_identity_transform()` may no longer be needed for this path, because BayesTools returns `NULL` for identity transforms.

Suggested tests:

- Standardized continuous moderator prior on original scale, base and ggplot.
- Standardized scale-regression prior on original scale, base and ggplot.
- Spike-and-slab/point-mass prior still plots.
- Identity-transform case returns to RoBMA's ordinary prior plotting fallback.

### Fitted formula design metadata

Use:

```r
BayesTools::JAGS_formula_design(fit, parameter = NULL)
```

Storage:

```r
attr(fit, "formula_design")[[parameter]]
```

Each design currently includes:

```r
list(
  parameter,
  formula,
  model_frame,
  model_matrix,
  column_names,
  raw_column_names,
  assign,
  terms,
  contrasts,
  xlevels,
  predictors,
  predictor_types,
  model_terms,
  model_terms_type,
  prior_list,
  formula_scale,
  rank,
  qr_pivot,
  aliased,
  transformed_terms,
  random_effects,
  jags_data_names
)
```

Important details:

- `model_matrix` is the exact design matrix BayesTools used to construct JAGS data after formula processing, predictor scaling, contrast assignment, and JAGS-safe column-name conversion.
- `column_names` are JAGS-safe names. Interactions use `__xXx__`.
- `raw_column_names` keep the original `model.matrix()` names before JAGS-safe conversion.
- `rank`, `qr_pivot`, and `aliased` are computed from the stored `model_matrix`.
- The accessor returns `NULL` if no design metadata exists.
- With `parameter = NULL`, it returns the full named list.
- With an unknown parameter, it errors.
- Design metadata is attached by `JAGS_fit()` and preserved by `JAGS_extend()`.
- The formula stored in `formula_design$formula` is metadata for the design matrix. Do not use its environment for evaluation. Use `model_frame`, `model_matrix`, and `terms` from the design object instead.

RoBMA connection:

- In `R/residuals.R::.get_model_matrix()`, replace moderator `model.matrix()` reconstruction with the BayesTools design for `"mu"` when available.
- Keep the local one-column intercept fallback for intercept-only models without a formula.
- Keep local PET/PEESE augmentation because those are RoBMA bias terms, not BayesTools formula terms.
- After appending PET/PEESE columns, recompute the final rank with `qr(X)$rank` for diagnostics using the augmented matrix. For the moderator-only part, prefer `design$rank`.

Suggested implementation shape:

```r
if (!is_mods) {
  K <- nrow(object[["data"]][["outcome"]])
  X <- matrix(1, nrow = K, ncol = 1)
  colnames(X) <- "(Intercept)"
  attr(X, "assign") <- 0L
} else {
  design <- BayesTools::JAGS_formula_design(object[["fit"]], "mu")

  if (!is.null(design)) {
    X <- design[["model_matrix"]]
  } else {
    mods_data    <- object[["data"]][["mods"]]
    mods_formula <- attr(mods_data, "formula")
    X <- model.matrix(mods_formula, data = mods_data)
  }
}
```

If user-facing diagnostics need display names, map `colnames(X)` through:

```r
stats::setNames(design[["raw_column_names"]], design[["column_names"]])
```

Suggested tests:

- Factor contrasts match the fitted BayesTools design.
- Standardized continuous moderator columns match BayesTools design values.
- Interactions use the JAGS-safe names internally.
- Rank-deficient moderator designs expose `design$rank < ncol(design$model_matrix)` and `any(design$aliased)`.
- Existing intercept-only and PET/PEESE diagnostics still work.

### Weightfunction observed p-value rug

Use:

```r
BayesTools::plot_posterior(
  samples,
  parameter = "weightfunction",
  ...,
  data = NULL,
  show_data = FALSE,
  dots_data = list()
)
```

Behavior:

- Supported only for full weightfunction posterior plots, not individual component plots.
- `data` can be a numeric p-value vector in `[0, 1]` or a data frame with a `p` column.
- Base and ggplot output both draw observed p-values as rug ticks.
- `rescale_x = TRUE` is respected by mapping p-values onto the same transformed x-axis as the weightfunction.
- `dots_data` supports at least:
  - `col` or `color`
  - `alpha`
  - `lwd`, `linewidth`, or `size`
  - `side` or `rug_side`
  - `height`, `rug_height`, or `ticksize`

RoBMA connection:

- In `R/plot.R::plot_weightfunction.brma()`, remove the current hard stop for `show_data`.
- Accept and validate `show_data` and optionally `dots_data`.
- Compute observed p-values with the same signed one-sided rule as selection mapping.

Suggested implementation shape:

```r
dots_raw <- list(...)
show_data <- isTRUE(dots_raw[["show_data"]])
dots_data <- dots_raw[["dots_data"]]
if (is.null(dots_data)) {
  dots_data <- list()
}
dots_raw[["show_data"]] <- NULL
dots_raw[["dots_data"]] <- NULL

BayesTools::check_bool(show_data, "show_data")
BayesTools::check_list(dots_data, "dots_data")

data_p <- NULL
if (show_data) {
  selection_context <- .selection_context(x)
  data_p <- stats::pnorm(
    selection_context[["sign"]] *
      selection_context[["yi"]] / selection_context[["sei"]],
    lower.tail = FALSE
  )
}

args$show_data <- show_data
args$data      <- data_p
args$dots_data <- dots_data
```

Then continue to call:

```r
plot <- suppressMessages(do.call(BayesTools::plot_posterior, args))
```

Suggested tests:

- `plot_weightfunction(fit, show_data = TRUE)` no longer errors.
- Base plot with `show_data = TRUE`.
- Ggplot with `show_data = TRUE`.
- `rescale_p_values = TRUE` and `FALSE`.
- Styling pass-through via `dots_data = list(color = "red", alpha = .5, linewidth = .4, rug_side = "top")`.

## Dependency/version handling

The local BayesTools `DESCRIPTION` still says `Version: 0.3.0`. Until BayesTools is version-bumped, RoBMA should use feature checks instead of only a version requirement.

Suggested guard:

```r
.check_bayestools_forward_api <- function() {
  required <- c(
    "formula_add_intercept",
    "plot_transformed_prior",
    "JAGS_formula_design"
  )

  missing <- required[!vapply(
    required,
    function(x) exists(x, envir = asNamespace("BayesTools"), inherits = FALSE),
    logical(1)
  )]

  if (length(missing) > 0) {
    stop(
      "RoBMA requires a BayesTools build with forward APIs: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
```

## BayesTools verification already run

In `C:\R-Packages\BayesTools`:

- Roxygen regenerated `NAMESPACE` and Rd files.
- Focused regression tests passed for:
  - `tests/testthat/test-JAGS-formula.R`
  - `tests/testthat/test-JAGS-fit-edge-cases.R`
  - `tests/testthat/test-JAGS-fit.R`
  - `tests/testthat/test-priors-linear-density.R`
  - `tests/testthat/test-model-averaging-plots-edge-cases.R`
  - `tests/testthat/test-weightfunction-plot-analytic.R`
- `devtools::check(document = FALSE, build_args = "--no-build-vignettes", args = c("--no-manual", "--ignore-vignettes"))` completed with:
  - `0 errors`
  - `0 warnings`
  - `3 notes`

The check notes were unrelated to this API work: existing long snapshot file paths, time verification, and existing top-level files.

## Summary for the RoBMA agent

Do these four rewires:

1. Replace `BayesTools:::.add_intercept_to_formula()` with `BayesTools::formula_add_intercept()`.
2. Replace private transformed-prior density/plot helper calls with `BayesTools::plot_transformed_prior()`.
3. Replace moderator model-matrix reconstruction with `BayesTools::JAGS_formula_design(object[["fit"]], "mu")$model_matrix` when present.
4. Pass observed p-values to `BayesTools::plot_posterior(..., parameter = "weightfunction", show_data = TRUE, data = p, dots_data = dots_data)`.

After those changes, RoBMA should no longer need the private BayesTools calls listed in `docs/bayestools-forward-spec.md`.
