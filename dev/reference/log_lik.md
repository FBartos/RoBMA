# Extract Pointwise Log-Likelihood Draws

Extract posterior pointwise log-likelihood draws from a fitted `brma`
object. The result is an \\S \times K\\ or \\S \times G\\ matrix, where
\\S\\ is the number of posterior draws, \\K\\ is the number of
estimates, and \\G\\ is the number of clusters.

## Usage

``` r
log_lik(object, ...)

# S3 method for class 'brma'
log_lik(object, unit = "estimate", ...)
```

## Arguments

- object:

  a brma model object.

- ...:

  currently unused.

- unit:

  output unit. See
  [`add_loo`](https://fbartos.github.io/RoBMA/dev/reference/add_loo.brma.md).

## Value

An \\S \times K\\ or \\S \times G\\ matrix of pointwise log-likelihood
draws with `RoBMA_target` metadata.

## Details

This function returns posterior log-likelihood draws, not a scalar
maximized log-likelihood. The
[`logLik.brma`](https://fbartos.github.io/RoBMA/dev/reference/logLik.brma.md)
method returns the same matrix for the
[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) generic;
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) are not defined for it.

The log-likelihood is computed for each target at each posterior sample.
For binomial and Poisson models, each observation consists of a pair of
counts (ai/ci or x1i/x2i) that together define a single effect size
estimate. Estimate-unit columns condition on the retained observations
in each declared dependency block. For Gaussian targets after
integrating their declared local Gaussian sources, including correlated
known-`V`
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
and specialized multilevel normal models, these are Schur conditional
scores \\p(y_i \mid y\_{-i}, \theta)\\; their row sum generally is not
the full joint likelihood. Likelihood powers and their normalizing
constants are retained inside the local-effect integral for weighted
Gaussian cluster models. The observation CDF used by residuals remains
unpowered. Numerically unresolved retained moments raise
`RoBMA_weighted_gaussian_unavailable`. Likelihood powers and
retained-sampling context follow the declared target. Independent
observations reduce to ordinary pointwise contributions. Cluster-unit
columns delete each complete cluster jointly. Sampling-conditioned
selection models retain the declared remaining sampling errors and
random context for estimate deletion, including ordinary component
branches and constant selection weights. Cluster deletion keeps its
distinct joint target. Earlier affected estimate scores must be
recomputed. Estimate-unit GLMM draws integrate the nuisance base rate or
log-rate over its prior: point priors and untruncated beta (base rate)
and normal (log-rate) priors adaptively, other priors on a refined
prior-CDF grid; a uniform(0, 1) base-rate prior is a Beta(1, 1) prior
and is integrated adaptively. The adaptive Gauss-Hermite rules raise the
quadrature order until successive orders agree within the tolerance;
this agreement is an empirical refinement check, not an error bound. The
grid is accepted only when its log-likelihood changes by less than
`1e-6` between refinements. Any nuisance prior without an adaptive route
(other than an untruncated beta base rate or an untruncated normal
log-rate; truncated ones take the grid too) can fail to reach this,
especially with few events or very large arms, and then stops with the
error class `RoBMA_glmm_grid_nonconvergence` described below; gamma,
lognormal, and inverse-gamma log-rate priors and lognormal base-rate
priors are examples. A beta base-rate prior or a normal log-rate prior
avoids it. The automatic computation after fitting (the
`autocompute.loo` and `autocompute.waic` options) warns instead and
returns the fit without LOO/WAIC. Numerical quadrature failures share
class `RoBMA_glmm_numerical_error`. AGHQ failures have class
`RoBMA_glmm_aghq_nonconvergence`; prior-CDF refinement failures retain
`RoBMA_glmm_grid_nonconvergence`. Their fields identify the outcome
family, observation, and available sample/refinement diagnostics. A beta
baserate prior whose shapes are too extreme for the prior-CDF grid to be
represented in double precision is unavailable and raises
`RoBMA_glmm_prior_unavailable`, also of class
`RoBMA_glmm_numerical_error`. Cluster-unit GLMM draws are unavailable
until certified nested adaptive quadrature is implemented.

## See also

[`loo.brma`](https://fbartos.github.io/RoBMA/dev/reference/loo.brma.md),
[`add_waic`](https://fbartos.github.io/RoBMA/dev/reference/add_waic.brma.md)
