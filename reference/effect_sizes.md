# Effect size transformations

Functions for transforming between different effect size measures.

## Usage

``` r
d2r(d)

d2z(d)

d2logOR(d)

d2OR(d)

r2d(r)

r2z(r)

r2logOR(r)

r2OR(r)

z2r(z)

z2d(z)

z2logOR(z)

z2OR(z)

logOR2r(logOR)

logOR2z(logOR)

logOR2d(logOR)

logOR2OR(logOR)

OR2r(OR)

OR2z(OR)

OR2logOR(OR)

OR2d(OR)
```

## Arguments

- d:

  Cohen's d.

- r:

  correlation coefficient.

- z:

  Fisher's z.

- logOR:

  log(odds ratios).

- OR:

  offs ratios.

## Details

All transformations are based on (Borenstein et al. 2011) . In case that
a direct transformation is not available, the transformations are
chained to provide the effect size of interest.

## References

Borenstein M, Hedges LV, Higgins JP, Rothstein HR (2011). *Introduction
to meta-analysis*. John Wiley & Sons.

## See also

[`standard_errors()`](reference/standard_errors.md),
[`sample_sizes()`](reference/sample_sizes.md)
