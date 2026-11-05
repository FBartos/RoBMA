# Diagnose Selection Conditioning Sensitivity

Describes how integrating the contexts retained by a fitted selection
model would reweight their pre-selection Gaussian distribution. This
compares reporting models; it is not a numerical approximation error
bound.

## Usage

``` r
selection_sensitivity_diagnostics(
  object,
  max_posterior_samples = 256L,
  latent_samples = 512L,
  seed = 1L,
  integration_control = NULL
)

add_selection_sensitivity_diagnostics(object, ...)
```

## Arguments

- object:

  A fitted normal-outcome model with a bound selection model.

- max_posterior_samples:

  Maximum number of posterior draws to evaluate. A positive integer of
  at least `10`, or `Inf` for all posterior draws. Draws are
  deterministically thinned over posterior row order.

- latent_samples:

  Number of fresh pre-selection retained-context draws per posterior
  draw. Must be one finite integer of at least `2`; defaults to `512`.

- seed:

  Integer seed for the diagnostic simulation. The calling R session's
  random-number state is restored before returning.

- integration_control:

  Optional settings from
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/dev/reference/set_selection_likelihood_control.md)
  for diagnostic event integration. `NULL` uses the fitted model's
  settings.

- ...:

  Simulation arguments passed to `selection_sensitivity_diagnostics()`.

## Value

A data frame with one row per connected source unit, identified by
`block`, comma-separated observation `rows`, and `n_estimates`. The
`ess_fraction`, `total_variation`, and `log_weight_iqr` summaries each
have `_q05`, `_median`, and `_q95` columns across all evaluated
posterior draws, including non-selection draws. The `ess_fraction_mcse`
and `total_variation_mcse` columns have `_median` and `_max` suffixes
for conditional Monte Carlo uncertainty. `integration_relative_mcse_max`
and `integration_relative_quadrature_error_max` are maxima over
constituent event evaluations, not an error bound for the multiplied
unit weight. `posterior_samples`, `total_posterior_samples`, and
`latent_samples` report the simulation budgets. The `settings` attribute
records `max_posterior_samples`, `latent_samples`, `seed`, and
`integration_control`. The `target` attribute records the comparison,
fitted and reference models, applicability, source roles, fitted and
reference sampling structures, publication groups, likelihood unit, row
blocks, and resolved event integration settings. The add function
returns the fitted object.

## Details

At each posterior draw the fitted model defines a conditional Gaussian
response law, retained contexts, and any best-rule publication groups.
The diagnostic evaluates each Gaussian selection unit's normalizer at
fresh retained contexts, preserving the fitted product or best weight
rule and any best-rule publication groups. Product-only models require
no grouping. Within a connected source unit, their product gives the
relative reweighting when the retained contexts are integrated before
selection normalization. Sampling-conditioned comparisons use the fitted
candidate factors, including singular declared support. Integration
diagnostics report maxima over the constituent event evaluations, rather
than an error bound for their product.

`ess_fraction` estimates `mean(W)^2 / mean(W^2)`; this is a
retained-context reweighting ESS fraction within each posterior
parameter state, not the posterior sample ESS. `total_variation`
estimates the distance between the pre-selection context law and its
normalized `W`-reweighted law. The fitted and reference conditioning
specifications, publication identities, and source roles are retained in
the `target` attribute. With no retained context this comparison is
inapplicable and its reweighting is identically constant.

The fitted weight rule is held fixed. These summaries do not compare
arbitrary weight rules, identify the correct reporting mechanism, or
reveal parameter regions absent from the fitted posterior.

Simulation budgets control diagnostic Monte Carlo precision separately
from the fitted Gaussian event integration controls. The largest
unit-level posterior-median total variation distance is silent up to
`0.05`, produces a message above `0.05` and below `0.10`, and a warning
from `0.10`. These notifications describe sensitivity rather than
statistical adequacy.

This optional comparison is computed only when requested; fitting does
not compute it. Calls without simulation arguments reuse an explicitly
stored result, and explicit settings reuse it only when they match. Use
`add_selection_sensitivity_diagnostics()` to compute and store a result.
Chain extension silently refreshes only stored diagnostics using their
stored simulation settings; label-only updates retain them. Failed
refreshes remain stored as unavailable without discarding the fitted
posterior. Call `selection_sensitivity_diagnostics()` to inspect the
comparison or its error. Sensitivity notifications are emitted only by
these explicit helper calls, not by fitting, extension, or ordinary
model printing.

The substantial-sensitivity warning has class
`RoBMA_selection_sensitivity_warning`, retaining `simpleWarning`,
`warning`, and `condition`. Its fields are `maximum_total_variation`,
`threshold` (`0.10`), and `comparison`
(`"retained_context_reweighting"`). Own nonfinite normalizer refusals
have class `RoBMA_selection_sensitivity_normalizer_unavailable` and
fields `reason`, `posterior_draw` (the evaluated posterior-row index),
`latent_draw`, `event`, and the observed `log_mass` values. Unsupported
retained covariance has class
`RoBMA_selection_sensitivity_covariance_unavailable` and
`reason = "retained_covariance_unavailable"`. Both refusals inherit from
`RoBMA_selection_sensitivity_error`, `simpleError`, `error`, and
`condition`. Downstream integration failures retain their original
conditions.

## See also

[`bselmodel.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.mv.md),
[`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/dev/reference/set_selection_likelihood_control.md)
