# True Effects

Computes the estimated true effects (theta) from a fitted model. This is
a separate S3 generic whose `brma` method delegates to
[`blup.brma`](https://fbartos.github.io/RoBMA/reference/blup.brma.md).

## Usage

``` r
true_effects(object, ...)
```

## Arguments

- object:

  a fitted model object

- ...:

  additional arguments passed to methods

## Value

Method-specific return value, typically a summary table or posterior
samples of BLUP or empirical-Bayes true-effect summaries.

## See also

[`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md),
[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
