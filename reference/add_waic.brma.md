# Add WAIC to brma Objects

Compute the Widely Applicable Information Criterion (WAIC) for brma
model objects and store the result in the object.

## Usage

``` r
# S3 method for class 'brma'
add_waic(object, unit = "estimate", ...)
```

## Arguments

- object:

  a brma model object.

- unit:

  output/deletion unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md);
  the same accepted values and multilevel constraint apply.

- ...:

  additional arguments passed to
  [`waic`](https://mc-stan.org/loo/reference/waic.html).

## Value

The brma object with the WAIC result stored in
`object[["waic"]][[unit]]`.

## Details

WAIC is an alternative to LOO-CV for estimating out-of-sample predictive
accuracy. Like LOO, it evaluates expected predictive performance for new
observations.

In most cases, LOO-PSIS (via
[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md))
is preferred over WAIC because it provides better estimates and includes
diagnostics (Pareto k values) that indicate when the approximation may
be unreliable.

For Gaussian multilevel and known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models, estimate-unit WAIC uses the same integrated conditional
log-likelihood matrix as estimate-unit LOO, \\p(y_i \mid y\_{-i},
\theta)\\. Estimate deletion retains the fitted grouping and dependency
structure. Weighted Gaussian cluster models integrate the powered row
likelihood, with its normalizing constant, over the retained-data
cluster-effect distribution. Earlier stored weighted estimate-unit WAIC
results require recomputation with `add_waic()`.

For [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models with known random-effect group covariance `R`, the known `R`
matrix and fitted random-effect covariance enter through BayesTools'
metadata-defined marginal \\ZGZ'\\ covariance. The same four-field
comparison policy applies across sampled, marginalized, and mixed
local-effect representations.

## See also

[`waic.brma`](https://fbartos.github.io/RoBMA/reference/waic.brma.md),
[`add_loo`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md),
[`waic`](https://mc-stan.org/loo/reference/waic.html)
