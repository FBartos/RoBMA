# Standard errors transformations

Functions for transforming between standard errors of different effect
size measures.

## Usage

``` r
se_d2se_logOR(se_d, logOR)

se_d2se_r(se_d, d)

se_r2se_d(se_r, r)

se_logOR2se_d(se_logOR, logOR)

se_d2se_z(se_d, d)

se_r2se_z(se_r, r)

se_r2se_logOR(se_r, r)

se_logOR2se_r(se_logOR, logOR)

se_logOR2se_z(se_logOR, logOR)

se_z2se_d(se_z, z)

se_z2se_r(se_z, z)

se_z2se_logOR(se_z, z)
```

## Arguments

- se_d:

  standard error of Cohen's d

- logOR:

  log(odds ratios)

- d:

  Cohen's d

- se_r:

  standard error of correlation coefficient

- r:

  correlation coefficient

- se_logOR:

  standard error of log(odds ratios)

- se_z:

  standard error of Fisher's z

- z:

  Fisher's z

## Details

Transformations for Cohen's d, Fisher's z, and log(OR) are based on
(Borenstein et al. 2011) . Calculations for correlation coefficient were
modified to make the standard error corresponding to the computed on
Fisher's z scale under the same sample size (in order to make all other
transformations consistent). In case that a direct transformation is not
available, the transformations are chained to provide the effect size of
interest.

It is important to keep in mind that the transformations are only
approximations to the true values. From our experience, `se_d2se_z`
works well for values of se(Cohen's d) \< 0.5. Do not forget that the
effect sizes are standardized and variance of Cohen's d = 1. Therefore,
a standard error of study cannot be larger unless the participants
provided negative information (of course, the variance is dependent on
the effect size as well, and, can therefore be larger).

When setting prior distributions, do NOT attempt to transform a standard
normal distribution on Cohen's d (mean = 0, sd = 1) to a normal
distribution on Fisher's z with mean 0 and sd = `se_d2se_z(0, 1)`. The
approximation does NOT work well in this range of values. Instead,
approximate the sd of distribution on Fisher's z using samples in this
way: `sd(d2z(rnorm(10000, 0, 1)))` or, specify the distribution on
Cohen's d directly.

## References

Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2011). *Introduction
to meta-analysis*. John Wiley & Sons.

## See also

[`effect_sizes()`](https://fbartos.github.io/RoBMA/reference/effect_sizes.md),
[`sample_sizes()`](https://fbartos.github.io/RoBMA/reference/sample_sizes.md)
