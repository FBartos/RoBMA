# Update a brma Fit

Extends an existing fitted `brma` object by additional MCMC samples,
updates study labels, and optionally recomputes cached fit-dependent
quantities.

## Usage

``` r
# S3 method for class 'brma'
update(
  object,
  formula. = NULL,
  ...,
  sample_extend = NULL,
  adapt = NULL,
  slab = NULL,
  autofit_control = NULL,
  convergence_checks = NULL,
  recompute = c("all", "drop"),
  parallel = NULL,
  cores = NULL,
  silent = NULL,
  evaluate = TRUE
)
```

## Arguments

- object:

  a fitted `brma` object.

- formula.:

  unsupported; included for compatibility with
  [`update`](https://rdrr.io/r/stats/update.html).

- ...:

  unsupported additional arguments.

- sample_extend:

  number of additional samples per chain, or `"auto"` to size this one
  extension from the convergence check of `object` against its stored
  autofit criteria, as `set_autofit_control(sample_extend = "auto")`
  sizes automatic extensions (at most twice the current samples per
  chain, but at least 1000; `max_time` does not shorten it, as no draws
  of the call have been timed). Defaults to `NULL`, which does not
  extend the chains; the `sample_extend` stored in the autofit control,
  including its default `"auto"`, is not used then.

- adapt:

  integer. Number of adaptation iterations run when JAGS rebuilds the
  sampler for this extension. Defaults to `NULL`, which uses the `adapt`
  of the initial fit. A supplied value applies to this extension only
  and does not replace the initial setting that later extensions without
  `adapt` use. Requires `sample_extend`.

- slab:

  optional character vector of study labels. Updating labels does not
  refit or extend the model.

- autofit_control:

  list of autofit control settings. Values are merged with the existing
  settings before extending.

- convergence_checks:

  list of convergence check settings. Values are merged with the
  existing settings and used to re-check the fit.

- recompute:

  whether cached `loo`, `waic`, and `marglik` values already stored in
  `object` should be recomputed after extension (`"all"`) or dropped
  with a warning (`"drop"`).

- parallel:

  logical. Whether to extend chains in parallel.

- cores:

  integer. Number of cores to use when `parallel = TRUE`.

- silent:

  logical. Whether to suppress JAGS output during extension.

- evaluate:

  unsupported; included for compatibility with
  [`update`](https://rdrr.io/r/stats/update.html).

## Value

The updated `brma` object.

## Details

Extending a fit preserves the fitted model and its posterior samples and
adds new samples. No additional burn-in is requested. The extension
rebuilds the compiled JAGS model from the stored chain states, which
resets the sampler tuning, so JAGS adapts again for `adapt` iterations
(the initial fit's setting unless supplied) before collecting those
samples. Warnings that adaptation was incomplete are kept with the fit's
warnings. Retained selection caches store normalizer calculations, not
JAGS model pointers or sampler tuning. Prior, data, and model-structure
updates are intentionally not supported by this method. For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
objects, the stored `V`, `R`, known-`V` backend, and marginalized
random-effect metadata are preserved because the fitted model structure
is reused.

Extensions continue the random-number generator state stored by JAGS.
They cannot be restarted from a new seed without changing the meaning of
chain continuation, so [`update()`](https://rdrr.io/r/stats/update.html)
has no `seed` argument. Extending a JAGS 4 fit requires refitting the
original model under JAGS 5.

Retained selection-normalizer caches are restored when compatible. Set
`RoBMA.options(selection.cache_retain = TRUE)` to retain the cache after
extension, or use
[`remove_selection_cache()`](https://fbartos.github.io/RoBMA/dev/reference/remove_selection_cache.md)
to remove saved cache entries.

Only explicitly stored selection-sensitivity diagnostics are refreshed
after chain extension, silently and using their stored simulation
settings, independently of `recompute`. Models without a stored
comparison do not compute one. Label-only updates retain the stored
diagnostic. Use
[`selection_sensitivity_diagnostics()`](https://fbartos.github.io/RoBMA/dev/reference/selection_sensitivity_diagnostics.md)
to inspect the refreshed result or error.

Named `NULL` elements in `autofit_control` and `convergence_checks`
explicitly disable or clear the corresponding nullable setting. Omitted
elements, a top-level `NULL`, and an empty
[`list()`](https://rdrr.io/r/base/list.html) inherit the setting stored
in `object`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- update(fit, sample_extend = 1000)
fit <- update(fit, sample_extend = "auto")
fit <- update(fit, slab = paste("Study", seq_len(nobs(fit))))
} # }
```
