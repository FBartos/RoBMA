# Weighted multivariate normal distribution

Density function for the weighted multivariate normal distribution with
`mean`, covariance matrix `sigma`, critical values `crit_x`, and weights
`omega`.

## Arguments

- x:

  quantiles.

- p:

  vector of probabilities.

- mean:

  mean

- sigma:

  covariance matrix.

- crit_x:

  vector of critical values defining steps.

- omega:

  vector of weights defining the probability of observing a t-statistics
  between each of the two steps.

- type:

  type of weight function (defaults to `"two.sided"`).

- log, log.p:

  logical; if `TRUE`, probabilities `p` are given as `log(p)`.

## Value

`.dwmnorm_fast` returns a density of the multivariate weighted normal
distribution.

## See also

[Normal](https://rdrr.io/r/stats/Normal.html),
[weighted_normal](reference/weighted_normal.md)
