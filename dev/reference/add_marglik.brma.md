# Add Marginal Likelihood to brma Objects

Compute the marginal likelihood of a brma model using bridge sampling
and store the result in the object.

## Usage

``` r
# S3 method for class 'brma'
add_marglik(
  object,
  parallel = NULL,
  cores = NULL,
  repetitions = 1L,
  method = c("normal", "warp3"),
  maxiter = 10000L,
  silent = TRUE,
  ...
)
```

## Arguments

- object:

  a brma model object.

- parallel:

  whether bridge-density evaluations run in parallel. The default,
  `NULL`, inherits the setting used to fit `object`.

- cores:

  number of bridge-sampling worker processes. The default, `NULL`,
  inherits the fitted core count and is capped by
  `RoBMA.get_option("max_cores")`. It is used only when
  `parallel = TRUE`.

- repetitions:

  number of independent bridge-sampling repetitions.

- method:

  bridge transformation; either `"normal"` or `"warp3"`.

- maxiter:

  maximum number of bridge iterations per repetition.

- silent:

  whether bridge-sampling progress is suppressed. The typed
  iteration-limit diagnostic remains visible even when `silent = TRUE`.

- ...:

  reserved for future use. Unused arguments produce a warning.

## Value

The brma object with the marginal likelihood result stored in
`object[["marglik"]]`.

## Details

The marginal likelihood is computed using the `bridgesampling` package
via the
[`BayesTools::JAGS_bridgesampling`](https://fbartos.github.io/BayesTools/reference/JAGS_bridgesampling.html)
wrapper. The result is stored in the object and can be extracted using
[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/dev/reference/bridge_sampler.brma.md).
Proposal draws use the seed stored with the fit. When that seed is
non-NULL, the caller's random-number state is restored after the call; a
NULL fit seed uses and advances the caller's random-number stream. A
finite bridge estimate remains available when observed iteration counts
exceed `maxiter`. The stored result's `success` records finite
availability, not certified convergence. Its optional
`iteration_limit_exceeded` is `TRUE` for observed excess, `FALSE` for
complete counts without excess, and `NA` or absent when telemetry is
unavailable. Counts may include upstream retries. Observed excess raises
`BayesTools_marglik_iteration_limit`, inheriting from
`BayesTools_marglik_diagnostic`, under either `silent` setting. Inspect
`print(object[["marglik"]])` and its repetition diagnostics.
`bridge_sampler(object)` returns the original upstream bridge object; it
does not promise the stored wrapper's optional telemetry field.

Product-space model-averaging objects (`BMA.norm`, `BMA.glmm`, `BMA.mv`,
`RoBMA`, and `RoBMA.mv`) do not expose a bridge-sampling marginal
likelihood; use predictive comparison methods such as
[`loo.brma`](https://fbartos.github.io/RoBMA/dev/reference/loo.brma.md)
instead. For a single model whose bridge parameters are all fixed by
point priors, the log marginal likelihood is evaluated exactly at those
fixed values; no bridge-sampling repetitions are required. For
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
known-`V` objects, bridge sampling evaluates the joint likelihood
corresponding to the fitted known-`V` backend, not the conditional
estimate-wise target used by LOO/WAIC diagnostics. Sampled Gaussian
location random-effect blocks are integrated exactly during bridge
evaluation as their draw-specific \\ZGZ'\\ covariance. Their SD,
allocation, and correlation parameters remain bridge coordinates; only
the standardized latent effects are removed. This is a
reparameterization of the same marginal-likelihood target, not a
likelihood approximation. Fitted estimate-level marginalized blocks
remain in the diagonal row variance. For Gaussian models fitted with the
specialized `cluster` interface, the standardized cluster effects are
likewise integrated exactly. The bridge retains total heterogeneity and
its allocation while evaluating the implied
diagonal-plus-cluster-rank-one covariance, including row-specific scale
regression and likelihood weights. Selection models retain every context
that the fitted model conditions on, including its normalized prior, and
integrate only the Gaussian sources declared by that model. The
likelihood uses the same full publication events and integration plan as
fitting. A context integrated before selection normalization cannot be
replaced by a retained-context parameterization. With
`known_sampling_variance = "condition"`, the bridge retains the
independent Gaussian augmentation coordinates and their normalized
priors. The complete sampling error and random effects are reconstructed
for each proposed bridge state before evaluating its selection
normalizer.

A non-finite bridge repetition raises a
`RoBMA_marglik_repetition_failure` error, inheriting from
`RoBMA_marglik_error`. Its `parent` retains the complete upstream
condition, including the failed `repetitions` table. Increase `maxiter`
when a failed repetition reached its iteration limit; otherwise inspect
the posterior samples and priors for non-finite densities. Fitting more
posterior draws (for example with a larger `sample` or
[`update()`](https://rdrr.io/r/stats/update.html)) may help in either
case. These actions do not guarantee that the diagnostics will pass.

A checked posterior sample that cannot form a full-rank bridge proposal
raises `RoBMA_marglik_rank_deficiency`, inheriting from
`RoBMA_marglik_error`. Its unchanged upstream `parent` and exposed
`reason`, `rank`, `n_coordinates`, `n_draws`, and
`dependent_coordinates` describe the observed sample. For
`reason = "insufficient_draws"`, fit more posterior draws with a larger
`sample` or [`update()`](https://rdrr.io/r/stats/update.html);
dependence may persist. For `reason = "dependent_coordinates"`, inspect
the parameterization and sampling; more draws alone need not resolve
dependence. The message names the first reported pivot coordinates; they
identify one dependent set, not uniquely faulty nodes.

The stored result records an identity of the data, weights, priors, and
model it was computed for.
[`bridge_sampler()`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html),
[`logml()`](https://rdrr.io/pkg/bridgesampling/man/logml.html),
[`bf()`](https://rdrr.io/pkg/bridgesampling/man/bf.html), and
[`post_prob()`](https://rdrr.io/pkg/bridgesampling/man/post_prob.html)
refuse a stored result that no longer matches the object, or that an
earlier RoBMA version stored, with an error of class
`RoBMA_marglik_stale` (inheriting from `RoBMA_marglik_error`); recompute
it with `add_marglik()`. `update(slab = )` keeps it, and extending the
fit with [`update()`](https://rdrr.io/r/stats/update.html) recomputes or
drops it.

## See also

[`bridge_sampler.brma`](https://fbartos.github.io/RoBMA/dev/reference/bridge_sampler.brma.md),
[`logml.brma`](https://fbartos.github.io/RoBMA/dev/reference/logml.brma.md),
[`bf.brma`](https://fbartos.github.io/RoBMA/dev/reference/bf.brma.md),
[`post_prob.brma`](https://fbartos.github.io/RoBMA/dev/reference/post_prob.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- brma(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit <- add_marglik(fit)

  bridge <- bridge_sampler(fit)
  print(bridge)

  logml(fit)
}
} # }
```
