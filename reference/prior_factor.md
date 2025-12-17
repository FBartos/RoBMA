# Creates a prior distribution for factors

`prior_factor` creates a prior distribution for fitting models with
factor predictors. (Note that results across different operating systems
might vary due to differences in JAGS numerical precision.)

## Usage

``` r
prior_factor(
  distribution,
  parameters,
  truncation = list(lower = -Inf, upper = Inf),
  prior_weights = 1,
  contrast = "meandif"
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

- contrast:

  type of contrast for the prior distribution. The possible options are

  `"meandif"`

  :   for contrast centered around the grand mean with equal marginal
      distributions, making the prior distribution exchangeable across
      factor levels. In contrast to `"orthonormal"`, the marginal
      distributions are identical regardless of the number of factor
      levels and the specified prior distribution corresponds to the
      difference from grand mean for each factor level. Only supports
      `distribution = "mnormal"` and `distribution = "mt"` which
      generates the corresponding multivariate normal/t distributions.

  `"orthonormal"`

  :   for contrast centered around the grand mean with equal marginal
      distributions, making the prior distribution exchangeable across
      factor levels. Only supports `distribution = "mnormal"` and
      `distribution = "mt"` which generates the corresponding
      multivariate normal/t distributions.

  `"treatment"`

  :   for contrasts using the first level as a comparison group and
      setting equal prior distribution on differences between the
      individual factor levels and the comparison level.

  `"independent"`

  :   for contrasts specifying dependent prior distribution for each
      factor level (note that this leads to an overparameterized model
      if the intercept is included).

## Value

return an object of class 'prior'.

## See also

[`prior()`](https://fbartos.github.io/BayesTools/reference/prior.html)

## Examples

``` r
# create an orthonormal prior distribution
p1 <- prior_factor(distribution = "mnormal", contrast = "orthonormal",
                   parameters = list(mean = 0, sd = 1))
```
