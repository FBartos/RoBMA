# Add LOO-PSIS to brma Objects

Compute approximate leave-one-out cross-validation (LOO-CV) using Pareto
smoothed importance sampling (PSIS) for brma model objects and store the
result in the object.

## Usage

``` r
# S3 method for class 'brma'
add_loo(object, unit = "estimate", r_eff = NULL, parallel = FALSE, ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output/deletion unit. `"estimate"` computes one contribution per
  effect-size estimate. `"cluster"` computes one contribution per
  cluster and is available only for multilevel models fitted through the
  specialized clustered interface.

- r_eff:

  optional vector of relative effective sample sizes. If not provided,
  it is computed from the log-likelihood values.

- parallel:

  Logical. If `TRUE`,
  [`loo::relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.html)
  and [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) use
  `RoBMA.get_option("max_cores")`. Log-likelihood construction is
  unchanged. If `FALSE`, those computations use one core.

- ...:

  additional arguments (currently ignored).

## Value

The brma object with the LOO result stored in `object[["loo"]][[unit]]`.

## Details

With `unit = "estimate"`, LOO-CV is computed with one contribution per
effect-size estimate. For binomial and Poisson models, each pair of
counts (ai/ci or x1i/x2i) that defines a single effect size estimate is
treated as one contribution. Estimate-unit GLMM log-likelihoods
integrate the nuisance base rate or log-rate over its prior: point
priors and untruncated beta (base rate) and normal (log-rate) priors
adaptively, other priors on a refined prior-CDF grid; a uniform(0, 1)
base-rate prior is a Beta(1, 1) prior and is integrated adaptively. The
grid is accepted only when its log-likelihood changes by less than
`1e-6` between refinements. Any nuisance prior without an adaptive route
(other than an untruncated beta base rate or an untruncated normal
log-rate; truncated ones take the grid too) can fail to reach this,
especially with few events or very large arms, and then stops with an
error of class `RoBMA_glmm_grid_nonconvergence` (also of class
`RoBMA_glmm_numerical_error`); gamma, lognormal, and inverse-gamma
log-rate priors and lognormal base-rate priors are examples. A beta
base-rate prior or a normal log-rate prior avoids it. When the
likelihood is computed automatically after fitting (the
`autocompute.loo` and `autocompute.waic` options), such a failure gives
a warning and the fit is returned without LOO/WAIC; `add_loo()` then
gives the same error.

With `unit = "cluster"`, LOO-CV is computed with one joint contribution
per cluster. For normal models without selection, weighted or not, this
uses the closed-form Gaussian integral over the held-out cluster effect,
with the likelihood weights applied as powers of the row densities.
Cluster-unit binomial and Poisson GLMM log-likelihoods are unavailable
until certified nested adaptive quadrature is implemented; use
`unit = "estimate"` for GLMMs. Cluster-unit contributions are also
unavailable for
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, whose random-formula deletion semantics still need a dedicated
design; use `unit = "estimate"` there.

Selection-model deletion conditions the Gaussian law on retained
outcomes and on the context declared by the fitted selection model. It
preserves the original publication event and selection standard errors.
Product weights on retained outcomes cancel; best selection continues to
depend on the best p-value among both retained and deleted outcomes.
This finite-vector deletion score differs from prediction for an
entirely new publication event.

For Gaussian multilevel and known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, estimate-unit LOO integrates Gaussian local effects and uses
\\p(y_i \mid y\_{-i}, \theta)\\ within each fitted dependency block.
Thus deleting one estimate retains the other estimates in its cluster or
random- effect block; it does not represent prediction for a wholly new
group. Gaussian `brma(..., cluster = ...)` models with non-unit
likelihood `weights` also integrate the cluster effect. Their score
integrates the entire powered row density against the cluster-effect
distribution given the retained weighted observations, including the
power-dependent normalizing constant. It is not the power of an
unweighted marginal density. Earlier stored weighted estimate-unit LOO
results must be recomputed with `add_loo()`; the fitted posterior does
not need refitting for this post-fit calculation. Selection and GLMM
models retain their supported routes. Numerically unresolved weighted
Gaussian retained moments stop with
`RoBMA_weighted_gaussian_unavailable`, a `RoBMA_numerical_unavailable`
condition.

For [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models with known random-effect group covariance `R`, the known `R`
matrix and fitted random-effect covariance enter the marginal \\ZGZ'\\
covariance through BayesTools' compiled random-effect metadata. The same
calculation applies to every supported random structure.

Estimate-unit deletion retains the fitted grouping and known-`V`
dependency structure; it does not represent deletion of an entire new
random group. Models with no local effects, sampled local effects,
marginalized local effects, or a mixture of sampled and marginalized
effects remain comparable when their data hash, deletion unit, retained
context, and likelihood target agree. The sampled/marginalized labels
are retained as provenance rather than comparison keys.

The PSIS object is essential for model comparison via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html) and
is automatically saved in the loo result. RoBMA stores target metadata
so comparisons can reject mismatched data, unit, or retained-context
targets. Each finite log-likelihood column observed to be constant
across posterior draws uses exact uniform importance ratios. Such
columns bypass relative-ESS and generalized-Pareto fitting, including
when other columns vary. Their Pareto-k diagnostic is recorded as the
sentinel zero, meaning constant ratios with no tail instability rather
than an estimated Pareto shape. Finite varying columns are centered by a
column-specific constant before relative-ESS and PSIS evaluation. The
offset is restored exactly afterward, preserving the LOO target while
avoiding probability-scale underflow.

**Important for model comparison:** When comparing models via
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html), the
selection is based on expected predictive performance for the documented
deletion target, not in-sample fit. For dependent estimates, deleting
one estimate may retain other observations from its fitted dependency
block; this is not new-group prediction.

## References

(Vehtari et al. 2017) (Vehtari et al. 2024)

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/reference/loo.brma.md),
[`loo`](https://mc-stan.org/loo/reference/loo.html),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html),
[`pareto_k_ids`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")
  fit <- bPET(yi = yi, vi = vi, data = dat.lehmann2018, measure = "SMD")

  fit <- add_loo(fit)
  loo_fit <- loo(fit)
  print(loo_fit)
  plot(loo_fit)
}
} # }
```
