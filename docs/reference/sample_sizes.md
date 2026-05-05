# Sample sizes to standard errors calculations

Functions for transforming between standard errors and sample sizes
(assuming equal sample sizes per group).

## Usage

``` r
se_d(d, n)

n_d(d, se)

se_r(r, n)

n_r(r, se)

se_z(n)

n_z(se)
```

## Arguments

- d:

  Cohen's d

- n:

  sample size of the corresponding effect size

- se:

  standard error of the corresponding effect size

- r:

  correlation coefficient

## Details

Calculations for Cohen's d, Fisher's z, and log(OR) are based on
(Borenstein et al. 2011) . Calculations for correlation coefficient were
modified to make the standard error corresponding to the computed on
Fisher's z scale under the same sample size (in order to make all other
transformations consistent). In case that a direct transformation is not
available, the transformations are chained to provide the effect size of
interest.

Note that sample size and standard error calculation for log(OR) is not
available. The standard error is highly dependent on the odds within the
groups and sample sizes for individual events are required.
Theoretically, the sample size could be obtained by transforming the
effect size and standard error to a different measure and obtaining the
sample size using corresponding function, however, it leads to a very
poor approximation and it is not recommended.

## References

Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2011). *Introduction
to meta-analysis*. John Wiley & Sons.

## See also

[`effect_sizes()`](https://fbartos.github.io/RoBMA/reference/effect_sizes.md),
[`standard_errors()`](https://fbartos.github.io/RoBMA/reference/standard_errors.md)
