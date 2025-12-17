# Creates a prior distribution for PET or PEESE models

`prior` creates a prior distribution for fitting a PET or PEESE style
models in RoBMA. The prior distribution can be visualized by the `plot`
function.

## Usage

``` r
prior_PEESE(
  distribution,
  parameters,
  truncation = list(lower = 0, upper = Inf),
  prior_weights = 1
)
```

## Arguments

- distribution:

  name of the prior distribution. The possible options are

  `"point"`

  :   for a point density characterized by a `location` parameter.

  `"normal"`

  :   for a normal distribution characterized by a `mean` and `sd`
      parameters.

  `"lognormal"`

  :   for a lognormal distribution characterized by a `meanlog` and
      `sdlog` parameters.

  `"cauchy"`

  :   for a Cauchy distribution characterized by a `location` and
      `scale` parameters. Internally converted into a generalized
      t-distribution with `df = 1`.

  `"t"`

  :   for a generalized t-distribution characterized by a `location`,
      `scale`, and `df` parameters.

  `"gamma"`

  :   for a gamma distribution characterized by either `shape` and
      `rate`, or `shape` and `scale` parameters. The later is internally
      converted to the `shape` and `rate` parametrization

  `"invgamma"`

  :   for an inverse-gamma distribution characterized by a `shape` and
      `scale` parameters. The JAGS part uses a 1/gamma distribution with
      a shape and rate parameter.

  `"beta"`

  :   for a beta distribution characterized by an `alpha` and `beta`
      parameters.

  `"exp"`

  :   for an exponential distribution characterized by either `rate` or
      `scale` parameter. The later is internally converted to `rate`.

  `"uniform"`

  :   for a uniform distribution defined on a range from `a` to `b`

- parameters:

  list of appropriate parameters for a given `distribution`.

- truncation:

  list with two elements, `lower` and `upper`, that define the lower and
  upper truncation of the distribution. Defaults to
  `list(lower = -Inf, upper = Inf)`. The truncation is automatically set
  to the bounds of the support.

- prior_weights:

  prior odds associated with a given distribution. The value is passed
  into the model fitting function, which creates models corresponding to
  all combinations of prior distributions for each of the model
  parameters and sets the model priors odds to the product of its prior
  distributions.

## Value

`prior_PET` and `prior_PEESE` return an object of class 'prior'.

## See also

[`plot.prior()`](https://fbartos.github.io/BayesTools/reference/plot.prior.html),
[`prior()`](https://fbartos.github.io/BayesTools/reference/prior.html)

## Examples

``` r
# create a half-Cauchy prior distribution
# (PET and PEESE specific functions automatically set lower truncation at 0)
p1 <- prior_PET(distribution = "Cauchy", parameters = list(location = 0, scale = 1))

plot(p1)
```
