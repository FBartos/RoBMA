# Estimated Marginal Means

S3 generic for estimated marginal means. The `brma` method works with
fitted moderator models and stores the BayesTools marginal-inference
object for [`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
marginal_means(object, ...)
```

## Arguments

- object:

  a fitted model object.

- ...:

  additional arguments passed to methods.

## Value

A method-specific estimated marginal means object.

## See also

[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`summary.brma()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md),
[`regplot()`](https://fbartos.github.io/RoBMA/reference/regplot.md)
