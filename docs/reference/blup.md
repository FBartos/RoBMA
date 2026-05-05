# Best Linear Unbiased Predictions (BLUPs)

Computes the estimated true effects (theta) from a fitted model. These
correspond to Best Linear Unbiased Predictions (BLUPs) or empirical
Bayes estimates.

## Usage

``` r
blup(object, ...)
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

[`predict.brma()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
