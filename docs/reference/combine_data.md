# Combines different effect sizes into a common metric

`combine_data` combines different effect sizes into a common measure
specified in `transformation`. Either a data.frame `data` with columns
named corresponding to the arguments or vectors with individual values
can be passed.

## Usage

``` r
combine_data(
  d = NULL,
  r = NULL,
  z = NULL,
  logOR = NULL,
  OR = NULL,
  t = NULL,
  y = NULL,
  se = NULL,
  v = NULL,
  n = NULL,
  lCI = NULL,
  uCI = NULL,
  study_names = NULL,
  cluster = NULL,
  weight = NULL,
  data = NULL,
  transformation = "fishers_z",
  return_all = FALSE,
  ...
)
```

## Arguments

- d:

  a vector of effect sizes measured as Cohen's d / Hedges' g
  (standardized mean differences)

- r:

  a vector of effect sizes measured as correlations

- z:

  a vector of effect sizes measured as Fisher's z

- logOR:

  a vector of effect sizes measured as log odds ratios

- OR:

  a vector of effect sizes measured as odds ratios

- t:

  a vector of t/z-statistics

- y:

  a vector of unspecified effect sizes (note that effect size
  transformations are unavailable with this type of input)

- se:

  a vector of standard errors of the effect sizes

- v:

  a vector of variances of the effect sizes

- n:

  a vector of overall sample sizes

- lCI:

  a vector of lower bounds of confidence intervals

- uCI:

  a vector of upper bounds of confidence intervals

- study_names:

  an optional argument with the names of the studies

- cluster:

  an optional argument specifying dependency between the studies (for
  using a multilevel model). Defaults to `NULL` for studies being
  independent.

- weight:

  specifies likelihood weights of the individual estimates. Notes that
  this is an untested experimental feature.

- data:

  a data frame with column names corresponding to the variable names
  used to supply data individually

- transformation:

  transformation to be applied to the supplied effect sizes before
  fitting the individual models. Defaults to `"fishers_z"`. We highly
  recommend using `"fishers_z"` transformation since it is the only
  variance stabilizing measure and does not bias PET and PEESE style
  models. The other options are `"cohens_d"`, correlation coefficient
  `"r"` and `"logOR"`. Supplying `"none"` will treat the effect sizes as
  unstandardized and refrain from any transformations.

- return_all:

  whether data frame containing all filled values should be returned.
  Defaults to `FALSE`

- ...:

  additional arguments.

## Value

`combine_data` returns a data.frame.

## Details

The aim of the function is to combine different, already calculated,
effect size measures. In order to obtain effect size measures from raw
values, e.g, mean differences, standard deviations, and sample sizes,
use [escalc](https://wviechtb.github.io/metafor/reference/escalc.html)
function.

The function checks the input values and in transforming the input into
a common effect size measure in the following fashion:

1.  obtains missing standard errors by squaring variances

2.  obtains missing standard errors from confidence intervals (after
    transformation to Fisher's z scale for `d` and `r`).

3.  obtains missing sample sizes (or standard errors for logOR) from
    t-statistics and effect sizes

4.  obtains missing standard errors from sample sizes and effect sizes

5.  obtains missing sample sizes from standard errors and effect sizes

6.  obtains missing t-statistics from sample sizes and effect sizes (or
    standard errors and effect sizes for logOR)

7.  changes the effect sizes direction to be positive

8.  transforms effect sizes into the common effect size

9.  transforms standard errors into the common metric

If the `transforms` is `NULL` or an unstandardized effect size `y` is
supplied, steps 4-9 are skipped.

## See also

[`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md),
[`check_setup()`](https://fbartos.github.io/RoBMA/reference/check_setup.md),
[`effect_sizes()`](https://fbartos.github.io/RoBMA/reference/effect_sizes.md),
[`standard_errors()`](https://fbartos.github.io/RoBMA/reference/standard_errors.md),
and
[`sample_sizes()`](https://fbartos.github.io/RoBMA/reference/sample_sizes.md)
