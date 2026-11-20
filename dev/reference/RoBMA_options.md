# Options for the RoBMA package

Inspect or change package-level defaults used by RoBMA.

## Usage

``` r
RoBMA.options(...)

RoBMA.get_option(name)
```

## Arguments

- ...:

  named option(s) to change. Names must be exact public option names,
  nonempty, and unique.

- name:

  a single non-missing character string matching one public option
  exactly; for available options, see details below.

## Value

`RoBMA.options()` invisibly returns a named list with all current
options after applying any changes. `RoBMA.get_option()` returns the
current value of the requested option.

## Details

The available options are:

- `max_cores`:

  number of cores to use for parallel computing (default is one fewer
  than detected logical cores, with a minimum/fallback of 1)

- `native_threads`:

  maximum number of native threads used to evaluate independent
  posterior rows of the compiled selection-likelihood batches in
  post-processing such as likelihood-aware densities, z-plots, and
  predictions. The default `NA` inherits the model's fitting setup:
  models fitted with `parallel = TRUE` evaluate their post-processing
  rows with up to `max_cores` native threads, unless this option sets an
  explicit thread count. Parallel work is arranged so that every
  reduction sums in a fixed order, so results do not depend on the
  thread count; single-state fitting and bridge-sampling evaluations
  stay serial. A model fitted with `parallel = FALSE` does its
  post-processing serially unless this option sets an explicit thread
  count.

- `check_scaling`:

  whether to check scaling of predictors (default `TRUE`)

- `silent`:

  whether to suppress JAGS output when fitting or extending models
  (default `FALSE`); the model-averaging constructors
  [`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md),
  [`BMA()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.md),
  [`BMA.norm()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.md),
  [`BMA.glmm()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.glmm.md),
  [`RoBMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md),
  and
  [`BMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/BMA.mv.md)
  default to `silent = TRUE` and ignore the option

- `jags.worker_output`:

  file path for parallel JAGS worker stdout and stderr when fitting or
  extending. The parent directory must exist; workers append to the same
  file and messages may interleave. The default empty string disables
  capture. Set this option inside any background job that fits the
  model. It does not change sampling or numerical integration settings.

- `autocompute.loo`:

  whether to automatically compute LOO (default `FALSE`); a GLMM whose
  nuisance-prior likelihood the quadrature cannot certify warns and the
  fit is returned without it

- `autocompute.waic`:

  whether to automatically compute WAIC (default `FALSE`); a GLMM whose
  nuisance-prior likelihood the quadrature cannot certify warns and the
  fit is returned without it

- `autocompute.marglik`:

  whether to automatically compute marginal likelihood (default `FALSE`)

- `cluster_likelihood.n_gamma`:

  number of Gauss-Hermite nodes formerly used for weighted normal
  cluster-unit log-likelihoods (default `15`). Normal cluster-unit
  log-likelihoods, weighted or not, are now integrated in closed form,
  so the option is still accepted but no longer affects them;
  selected-normal likelihoods use their own quadrature schedules and
  integration controls.

- `selection.cache_max_bytes`:

  one shared storage budget for exact and coarse selection-normalizer
  caches across all chains in a fit. The default `"auto"` uses the
  smaller of 4 GiB and one quarter of currently available system RAM,
  resolved before fitting. Supply a nonnegative whole number of bytes to
  set an explicit total; `0` disables caching. Storage grows as needed.
  Chains receive equal shares, pooled within each worker process;
  sequential chains share the total.

- `selection.cache_retain`:

  whether to save selection-normalizer cache entries with the fitted
  object after fitting or extending (default `FALSE`). Retained entries
  survive [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and are
  restored when extending a compatible fit. See
  [`remove_selection_cache`](https://fbartos.github.io/RoBMA/dev/reference/remove_selection_cache.md)
  for memory costs and removal.

- `selection.sampler`:

  the samplers of selection and model-averaged fits, in three tiers.
  `"slice.coarse"` (the default) adds the corrected selection sampler:
  coarse slice proposals with a full-target Metropolis correction for
  eligible scalar updates involving multivariate, cluster or factor step
  selection blocks (see
  [`selection_sampler_info`](https://fbartos.github.io/RoBMA/dev/reference/selection_sampler_info.md)
  for the blocks that take part, the fallback to the exact slice sampler
  and the per-fit record); a block that `sampler.blocks` declares is
  updated jointly by the corrected block sampler when it has such a
  selection block among its children, and by the factor slice sampler of
  the BayesTools JAGS module otherwise; parameters that `sampler.gates`
  declares keep the gated sampler. `"slice"` uses the samplers JAGS
  assigns together with the gated and factor slice samplers of the
  BayesTools JAGS module where `sampler.gates` and `sampler.blocks`
  declare them. `"jags.default"` uses only the samplers JAGS assigns, in
  every model and whatever `sampler.gates` and `sampler.blocks` say (as
  `"default"` did in RoBMA 4.1.1). The posterior target is unchanged.
  Settings are captured when each model is compiled, and
  [`update()`](https://rdrr.io/r/stats/update.html) keeps those of the
  fit. The deprecated values `"coarse_corrected"` (now `"slice.coarse"`)
  and `"default"` (now `"jags.default"`) are accepted with a warning of
  class `RoBMA_deprecated_option_value` (family `RoBMA_deprecated`).

- `selection.coarse_min_acceptance`:

  acceptance threshold of the corrected sampler (default `0.90`, a
  number between `0` and `1`). Each corrected sampler whose acceptance
  rate during adaptation is below it, that made fewer than 50 proposals
  during adaptation, or whose surrogate never differed from the exact
  density during adaptation updates its parameter or block with the
  exact sampler for the rest of the fit, and a fit with at least 1,000
  counted proposals whose post-adaptation acceptance rate, or that of
  its worst sampler, is below it warns. A block of d parameters is
  judged by its acceptance rate to the power 1/d (see
  [`selection_sampler_info`](https://fbartos.github.io/RoBMA/dev/reference/selection_sampler_info.md)).
  `0` disables the fallback and the warning. Captured when each model is
  compiled, like the other sampler settings.

- `selection.coarse_grid`:

  named nonnegative grid steps
  `c(mean = .05, variance = .05, log_weight = .01)` of the dense
  multivariate step blocks; cluster and factor blocks evaluate their
  coarse rule at the actual parameters and ignore the grid. Mean spacing
  is its step times the largest original sampling SE in the block;
  covariance-diagonal spacing is its step times that SE squared. Each
  diagonal is rounded upward on a grid whose origin is its own squared
  original sampling SE, retaining all off-diagonal entries. Positive
  weights are rounded on the natural-log scale; zero weights remain
  zero. A zero step retains that coordinate exactly but does not disable
  the coarse rule budget or fallback. These settings affect proposals,
  not the corrected target.

- `selection.coarse_max_rules`:

  coarse quadrature rule budget (default `3`, minimum `3`). Dense
  multivariate blocks try at most this many rules from the beginning of
  the fitted quadrature schedule for a coarse anchor; unsupported or
  unsuccessful coarse integration uses the prescribed unit normalizer,
  without QMC fallback. Cluster and factor blocks use the one rule of
  this number in their own schedule (the last rule of a shorter
  schedule) at the actual parameters, and the exact density where that
  rule is unavailable. The full-target evaluator retains its ordinary
  numerical controls and diagnostics.

- `sampler.gates`:

  whether the parameters of the publication-bias branches of
  model-averaged fits
  ([`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md),
  [`RoBMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md))
  are updated by the gated slice sampler of the BayesTools JAGS module
  (default `TRUE`): while the bias indicator does not select a branch,
  its parameters are drawn directly from their prior instead of by slice
  sampling. `selection.sampler = "jags.default"` turns the gated sampler
  off.

- `sampler.blocks`:

  whether the global parameters of selection models with estimated
  heterogeneity are updated jointly (default `TRUE`), slice-sampling
  along the principal directions of their posterior covariance learned
  during adaptation: under `selection.sampler = "slice.coarse"` by the
  corrected block sampler when the block has a dense, cluster or factor
  selection block that takes part in the corrected sampler
  ([`selection_sampler_info`](https://fbartos.github.io/RoBMA/dev/reference/selection_sampler_info.md)),
  otherwise and under `"slice"` by the factor slice sampler of the
  BayesTools JAGS module. Single models
  ([`bselmodel()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.md),
  [`bselmodel.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.mv.md))
  update the effect, moderator, heterogeneity and selection-weight
  parameters together, the directly sampled scalar random-effect SDs
  (the SD of a random intercept, or the one SD shared by the levels of
  `id()`, `diag(hom = TRUE)`, `cs()`, `ar1()` and `car()`) and the
  scalar correlations of `cs()`, `hcs()`, `ar1()`, `har()` and `car()`
  included. Never block members: random-effect SD priors over several
  levels and unstructured correlation matrices, vector (mean-difference
  and orthonormal) factor terms, the unit-level latent parameters and
  random-effect coordinates, and every mixture, spike-and-slab, ordered
  and point prior. Model-averaged fits
  ([`RoBMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md))
  update only the allocation SD and weights of the heterogeneity
  allocation together. `selection.sampler = "jags.default"` turns both
  block samplers off. Neither option changes the model or its posterior
  distribution, only the draws; settings are captured when each model is
  fitted, and [`update()`](https://rdrr.io/r/stats/update.html) keeps
  those of the fit.

- `default_UISD.effect`:

  default scaling of the unit information standard deviation for the
  effect size parameter (default `0.5`)

- `default_UISD.heterogeneity`:

  default scaling of the unit information standard deviation for the
  heterogeneity parameter (default `0.25`)

- `default_UISD.mods`:

  default scaling of the unit information standard deviation for the
  moderators (default `0.25`)

- `default_UISD.scale`:

  default scaling of the unit information standard deviation for the
  scale parameter (default `0.5`)

- `default_lograte.sd`:

  default standard deviation of the estimate-specific midpoint log-rate
  prior in Poisson GLMMs (default `1`)

- `default_informed_priors.mods`:

  default scaling of informed priors for moderators (default `0.5`)

- `default_informed_priors.scale`:

  default scaling of informed priors for the scale parameter (default
  `0.5`)

- `default_bias_weightfunction.alpha`:

  default alpha for the weightfunction (default `1`)

- `default_bias_PET.scale`:

  default scale for the PET (default `1`)

- `default_bias_PEESE.scale`:

  default scale for the PEESE (default `5`)

Boolean options require scalar `TRUE` or `FALSE` values. `max_cores`
must be integer-like and at least 1; `cluster_likelihood.n_gamma` must
be integer-like and at least 3. `native_threads` accepts `NA` or a
positive integer. `selection.coarse_min_acceptance` must be a number
between 0 and 1. Scale and alpha defaults must be finite positive
numbers.
