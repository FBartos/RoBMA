# Creates a prior distribution for a weight function

`prior_weightfunction` creates a prior distribution for fitting a RoBMA
selection model. The prior can be visualized by the `plot` function.

## Usage

``` r
prior_weightfunction(distribution, parameters, prior_weights = 1)
```

## Arguments

- distribution:

  name of the prior distribution. The possible options are

  `"two.sided"`

  :   for a two-sided weight function characterized by a vector `steps`
      and vector `alpha` parameters. The `alpha` parameter determines an
      alpha parameter of Dirichlet distribution which cumulative sum is
      used for the weights omega.

  `"one.sided"`

  :   for a one-sided weight function characterized by either a vector
      `steps` and vector `alpha` parameter, leading to a monotonic
      one-sided function, or by a vector `steps`, vector `alpha1`, and
      vector `alpha2` parameters leading non-monotonic one-sided weight
      function. The `alpha` / `alpha1` and `alpha2` parameters determine
      an alpha parameter of Dirichlet distribution which cumulative sum
      is used for the weights omega.

- parameters:

  list of appropriate parameters for a given `distribution`.

- prior_weights:

  prior odds associated with a given distribution. The model fitting
  function usually creates models corresponding to all combinations of
  prior distributions for each of the model parameters, and sets the
  model priors odds to the product of its prior distributions.

## Value

`prior_weightfunction` returns an object of class 'prior'.

## Details

Constrained cases of weight functions can be specified by adding
".fixed" after the distribution name, i.e., `"two.sided.fixed"` and
`"one.sided.fixed"`. In these cases, the functions are specified using
`steps` and `omega` parameters, where the `omega` parameter is a vector
of weights that corresponds to the relative publication probability
(i.e., no parameters are estimated).

## See also

[`plot.prior()`](https://fbartos.github.io/BayesTools/reference/plot.prior.html)

## Examples

``` r
p1 <- prior_weightfunction("one-sided", parameters = list(steps = c(.05, .10), alpha = c(1, 1, 1)))

# the prior distribution can be visualized using the plot function
# (see ?plot.prior for all options)
plot(p1)
```
