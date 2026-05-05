# BayesTools handoff spec

Scope: this records remaining BayesTools needs and the local RoBMA bridge-sampling fix that replaced previous private `bridgesampling:::` calls.

## 1. `BayesTools:::.add_intercept_to_formula`

Current dependency:

- `R/fit.R::.create_fit_formula_list()` calls `BayesTools:::.add_intercept_to_formula()` when a scale formula removes the intercept.
- Reason: scale regression uses an exponential/log-intercept parameterization; RoBMA warns and forces the intercept back because scale models cannot omit it.

Requested BayesTools API:

```r
BayesTools::formula_add_intercept(formula)
```

Requirements:

- Exported.
- Accepts one formula object, returns a formula object.
- Preserves formula environment.
- Removes `- 1`, `+ 0`, and equivalent no-intercept encodings robustly.
- Keeps response-free and response-containing formulas valid.
- Does not alter other terms, interactions, specials, offsets, or transformed calls.

RoBMA-side change after API lands:

- Replace `BayesTools:::.add_intercept_to_formula(formula)` with `BayesTools::formula_add_intercept(formula)`.
- Add a small regression test for scale formula `~ x - 1`, `~ 0 + x`, and environment preservation.

## 2. Prior plotting helpers

Current dependency:

- `R/plot.R::.plot_prior_unstandardized()` uses private BayesTools helpers:
  - `BayesTools:::.generate_transformed_prior_densities()`
  - `BayesTools:::.prior_linear_density_to_plot_data()`
  - `BayesTools:::.plot_prior_list.both()`
- Reason: RoBMA standardizes continuous predictors before fitting, but users need prior plots on the original predictor scale.

Requested BayesTools API, preferred high-level wrapper:

```r
BayesTools::plot_transformed_prior_linear(
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
Update: make this a more general function, i.e., beyond linear transform. examine the options and make sure it is consistent with the rest of Bayestools API

Requirements:

- Uses the same density transformation as current private helpers.
- Supports continuous density and point-mass/spike priors in one plot.
- Returns invisible base-plot metadata for `plot_type = "base"` and a ggplot for `plot_type = "ggplot"`.
- Accepts `formula_scale` in the same shape as `attr(fit, "formula_scale")`.
- Returns `NULL` or a clear classed condition when the requested parameter is an identity transform, so RoBMA can fall back to regular `plot_prior()`.

Minimum fallback if a wrapper is too large:

- Export stable equivalents of the three private helpers above with documented inputs, outputs, and classes.

RoBMA-side change after API lands:

- Replace the three `BayesTools:::` calls in `.plot_prior_unstandardized()`.
- Keep RoBMA responsible for selecting `mods` vs `scale` priors and user-facing parameter labels.
- Add tests that standardized continuous moderator/scale priors plot on original scale for base and ggplot.

## 3. Fitted formula design / model matrix metadata

Current issue:

- `R/residuals.R::.get_model_matrix()` reconstructs `model.matrix()` from stored RoBMA data.
- The fitted design is actually created inside `BayesTools::JAGS_formula()` during `BayesTools::JAGS_fit()`.
- Rebuilding later can diverge for factor contrasts, standardized predictors, transformed/interacted terms, aliases, rank deficiency, and any future BayesTools formula processing.

Requested BayesTools API:

1. `BayesTools::JAGS_fit()` should store the fitted formula design for each formula parameter:

```r
attr(fit, "formula_design")[[parameter]]
```

2. Export an accessor:

```r
BayesTools::JAGS_formula_design(fit, parameter = NULL)
```

Return shape for each `parameter` (`"mu"`, `"log_tau"`, etc.):

```r
list(
  parameter       = "mu",
  formula         = <processed formula used by JAGS_formula>,
  model_frame     = <fitted model.frame>,
  model_matrix    = <exact fitted model.matrix after contrast/scaling processing>,
  column_names    = <JAGS-safe coefficient column names>,
  raw_column_names = <pre-sanitized model.matrix names, if different>,
  assign          = attr(model_matrix, "assign"),
  terms           = terms_object,
  contrasts       = attr(model_matrix, "contrasts"),
  xlevels         = <factor levels used during fit>,
  predictors      = <predictor names>,
  predictor_types = <continuous/factor>,
  model_terms     = <term labels, including intercept when present>,
  model_terms_type = <continuous/factor by term>,
  prior_list      = <post-processed formula prior list for this parameter>,
  formula_scale   = <scale info for this parameter, or NULL>,
  rank            = qr(model_matrix)$rank,
  qr_pivot        = qr(model_matrix)$pivot,
  aliased         = <named logical over columns>,
  transformed_terms = <parsed expression/transformation metadata, if any>,
  random_effects  = <parsed random-effect metadata, if any>,
  jags_data_names = <mapping of terms to generated JAGS data nodes>
)
```

Requirements:

- `model_matrix` must be the exact matrix used to create JAGS data, after predictor scaling and contrast assignment.
- Preserve JAGS-safe names (`:` -> `__xXx__`) and raw display names.
- Include enough rank/alias metadata for diagnostics to choose effective fixed-effect rank instead of `ncol(X)`.
- Include transformed terms and expression terms even if they do not map 1:1 to `model.matrix()` columns.
- For factor priors, expose the final factor levels, contrasts, cell names, and any prior attributes BayesTools uses for evaluation.
- Accessor must work on fitted objects and on failed `BayesTools_fit` objects when `JAGS_formula()` completed before sampling failure.

RoBMA-side change after API lands:

- Replace `.get_model_matrix()` reconstruction with `BayesTools::JAGS_formula_design(object[["fit"]], "mu")$model_matrix`.
- For intercept-only models without a formula, keep local one-column intercept fallback.
- Keep local PET/PEESE augmentation, because those are RoBMA bias terms, not BayesTools formula terms.
- Use `rank` from BayesTools design for Cook's distance / leverage degrees of freedom; warn or drop aliased columns when needed.
- Update fitted-design tests to compare contrasts, standardized continuous predictors, interactions, and aliased designs.

## 4. Weightfunction plot `show_data`

Current issue:

- `plot_weightfunction.brma(show_data = TRUE)` validates `show_data` but never passes data to BayesTools or draws it locally.
- RoBMA currently delegates the full weightfunction plot to `BayesTools::plot_posterior(parameter = "weightfunction")`.

Requested BayesTools API:

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

Requirements:

- When `parameter = "weightfunction"` and `show_data = TRUE`, draw observed p-values as rug/tick marks on the weightfunction x-axis.
- Support both base and ggplot output.
- Respect `rescale_x`; data p-values must be transformed by the same x-axis mapping as the weightfunction.
- Accept `data` as a numeric vector of p-values in `[0, 1]` or a data frame with at least `p`.
- `dots_data` should support at least color, alpha, line width, and rug side/height.

RoBMA-side change after API lands:

- Compute observed p-values from the same signed one-sided rule used in `.selection_spec()`:

```r
p <- stats::pnorm(selection_spec$sign * yi / sei, lower.tail = FALSE)
```

- Pass `data = p`, `show_data = show_data`, and `dots_data` through to `BayesTools::plot_posterior()`.
- Add base/ggplot visual tests and an `as_data`/plot-data test if BayesTools exposes plot data.

